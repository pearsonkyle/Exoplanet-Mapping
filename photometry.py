import os
import glob
import json
import pickle
import argparse
import urllib.request
from io import BytesIO
from collections import Counter
import pyvo as vo
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import matplotlib.pyplot as plt
from scipy.stats import mode
from scipy.optimize import least_squares
from photutils import aperture_photometry, CircularAperture, CircularAnnulus
from scipy.ndimage import binary_dilation

target_coords = {
    'hd209458':[330.79488, 18.8843, 29.766, -17.976], # ra, dec in deg, proper motion in mas/year for ra and dec
    'hd189733':[300.1823, 22.7108, -3.208, -250.323]
}

def parse_args():
    parser = argparse.ArgumentParser()

    help_ = "Choose a target to process (default = all)"
    parser.add_argument("-t", "--target", help=help_, type=str, default="all")

    parser.add_argument("-f", "--filter", 
        help="Choose a filter (3.6, 4.5, 5.8, 8.0)", type=str, default="all")

    parser.add_argument("-m", "--mode", 
        help="Choose a mode to process (all/sub/full)", type=str, default="all")

    parser.add_argument("-r", "--reset", 
        help="Reset file list in order to reprocess all data", action='store_true', default=False)

    help_ = "Directory containing data (default = DATA/)"
    parser.add_argument("-d", "--datadir", help=help_, default="DATA/", type=str)

    help_ = "AOR observation filter"
    parser.add_argument("-a", "--aor", help=help_, default=0, type=int)

    return parser.parse_args()

def mesh_box(pos,box):
    pos = [int(np.round(pos[0])),int(np.round(pos[1]))]
    x = np.arange(pos[0]-box, pos[0]+box+1)
    y = np.arange(pos[1]-box, pos[1]+box+1)
    xv, yv = np.meshgrid(x, y)
    return xv.astype(int),yv.astype(int)

def mixed_psf(x,y,x0,y0,a,sigx,sigy,rot,b, w):
    gaus = gaussian_psf(x,y,x0,y0,a,sigx,sigy,rot, 0)
    lore = lorentz_psf(x,y,x0,y0,a,sigx,sigy,rot, 0)
    return (1-w)*gaus + w*lore + b

def gaussian_psf(x,y,x0,y0,a,sigx,sigy,rot, b):
    rx = (x-x0)*np.cos(rot) - (y-y0)*np.sin(rot)
    ry = (x-x0)*np.sin(rot) + (y-y0)*np.cos(rot)
    gausx = np.exp(-(rx)**2 / (2*sigx**2))
    gausy = np.exp(-(ry)**2 / (2*sigy**2))
    return a*gausx*gausy + b

def lorentz_psf(x,y,x0,y0,a,sigx,sigy,rot, b):
    rx = (x-x0)*np.cos(rot) - (y-y0)*np.sin(rot)
    ry = (x-x0)*np.sin(rot) + (y-y0)*np.cos(rot)
    lorex = sigx**2 / (rx**2 + sigx**2)
    lorey = sigy**2 / (ry**2 + sigy**2)
    return a*lorex*lorey + b

def fit_psf(data, pos, init, lo, up, psf_function=gaussian_psf, lossfn='linear', method='trf', box=15):
    xv, yv = mesh_box(pos, box)

    def fcn2min(pars):
        model = psf_function(xv, yv, *pars)
        return (data[yv, xv] - model).flatten()

    if method == 'trf':
        res = least_squares(fcn2min, x0=[*pos, *init], bounds=[lo, up], loss=lossfn, jac='3-point', method='dogbox',
                            xtol=None, ftol=1e-3, tr_options='exact')
    else:
        res = least_squares(fcn2min, x0=[*pos, *init], loss=lossfn, jac='3-point', method=method)
    return res.x

def psf_phot(data, xc, yc, box=10):

    # only fit ROI
    xv, yv = mesh_box([xc,yc], box)

    # estimate psf shape
    init = [np.nanmax(data[yv, xv]) - np.nanmin(data[yv, xv]), 1, 1, 0, np.nanmin(data[yv, xv]), 0.5]

    pars = fit_psf(
        data,
        [xc, yc],  # position estimate
        init,  # initial guess: [amp, sigx, sigy, rotation, bg] 
        # bounds: [xc, yc, amp, sigx, sigy, rotation,  bg]
        [xc - box * 0.5, yc - box * 0.5, 0, 0.01, 0.01, -np.pi / 4, np.nanmin(data) - 1, 0], # lower bound
        [xc + box * 0.5, yc + box * 0.5, 1e7, 20, 20, np.pi / 4, np.nanmax(data[yv, xv]) + 1, 1],  # upper bound
        psf_function=mixed_psf, method='trf',
        box=box  # only fit a subregion +/- 5 px from centroid
    )

    #model = mixed_psf(xv, yv, *pars)
    #res = data[yv, xv] - model
    return pars # [xc, yc, amp, sigx, sigy, rotation,  bg]

def find_target(target, hdu, stats=[], verbose=False):

    if len(stats) == 0:
        # query simbad to get proper motions
        service = vo.dal.TAPService("http://simbad.u-strasbg.fr/simbad/sim-tap")
        # http://simbad.u-strasbg.fr/simbad/tap/tapsearch.html
        query = '''
        SELECT basic.OID, ra, dec, main_id, pmra, pmdec
        FROM basic JOIN ident ON oidref = oid
        WHERE id = '{}';
        '''.format(target)
        result = service.search(query)
        ra = result['ra'][0]
        dec = result['dec'][0]
        pmra = result['pmra'][0]
        pmdec = result['pmdec'][0]
    else:
        ra,dec,pmra,pmdec = stats

    # set up astropy object
    coord = SkyCoord(
        ra = ra*u.deg,
        dec = dec*u.deg,
        distance=1*u.pc, 
        pm_ra_cosdec=pmra*u.mas/u.yr,
        pm_dec=pmdec*u.mas/u.yr,
        frame="icrs",
        obstime=Time("2000-1-1T00:00:00")
    )

    # apply proper motion
    t = Time(hdu.header['DATE_OBS'], format='isot', scale='utc')
    coordpm = coord.apply_space_motion(new_obstime=t)

    # wcs coordinate translation
    try:
        wcs = WCS(hdu.header)
    except ValueError:
        # FIXME https://github.com/astropy/astropy/issues/10527
        hdu.header['NAXIS'] = 2
        wcs = WCS(hdu.header)

    pixcoord = wcs.wcs_world2pix([[coordpm.ra.value, coordpm.dec.value]],0)

    if verbose:
        print("Simbad:",result)
        print("\nObs Date:",hdu.header['DATE_OBS'])
        print("NEW:", coordpm.ra, coordpm.dec)
        print("Pixels:",pixcoord[0])

    return pixcoord[0]

def phot(data,xc,yc,r=5,dr=5):
    # aperture photometry with background subtraction
    if dr>0:
        bgflux = skybg_phot(data,xc,yc,r+1,dr)
    else:
        bgflux = 0
    positions = [(xc, yc)]
    data = data-bgflux
    data[data<0] = 0 

    apertures = CircularAperture(positions, r=r)
    phot_table = aperture_photometry(data, apertures, method='exact')
    try:
        subdata = data[
            max(0,apertures.to_mask()[0].bbox.iymin):min(apertures.to_mask()[0].bbox.iymax,data.shape[0]-1),
            max(0,apertures.to_mask()[0].bbox.ixmin):min(apertures.to_mask()[0].bbox.ixmax,data.shape[1]-1)
        ] * apertures.to_mask()[0].data
    except:
        # happens when star is near edge of image, sizes are clipped
        subdata = data[
            max(0,apertures.to_mask()[0].bbox.iymin):min(apertures.to_mask()[0].bbox.iymax,data.shape[0]-1),
            max(0,apertures.to_mask()[0].bbox.ixmin):min(apertures.to_mask()[0].bbox.ixmax,data.shape[1]-1)
        ] 
    npp = subdata.sum()**2 / np.sum(subdata**2)
    return float(phot_table['aperture_sum']), bgflux, npp, apertures.to_mask()[0].data.sum()

def skybg_phot(data,xc,yc,r=10,dr=5,ptol=95):    
    # create a crude annulus to mask out bright background pixels 
    xv,yv = mesh_box([xc,yc], np.round(r+dr) )
    xv = np.clip(xv,0,data.shape[1]-1)
    yv = np.clip(yv,0,data.shape[0]-1)
    rv = ((xv-xc)**2 + (yv-yc)**2)**0.5
    mask = (rv>r) & (rv<(r+dr))
    cutoff = np.nanpercentile(data[yv,xv][mask], ptol)
    dat = np.array(data,dtype=float)
    dat[dat>cutoff] = np.nan # ignore bright pixels like stars 
    return np.nanmedian(dat)


def calc_moments(data: np.ndarray) -> tuple:
    """Returns (height, x, y, width_x, width_y, skew_x, skew_y, kurt_x, kurt_y)
    the gaussian parameters of a 2D distribution by calculating its moments"""
    
    # Check that the input data is a 2D numpy ndarray
    if not isinstance(data, np.ndarray):
        raise ValueError("data should be a numpy ndarray")
    if data.size == 0:
        raise ValueError("data should not be empty")
    if data.ndim != 2:
        raise ValueError("data should be a 2D numpy ndarray")
    
    # Calculate the total value of the data
    total = data.sum()
    
    # Calculate the maximum value of the data
    height = data.max()
    
    # Calculate the indices of the data
    X, Y = np.indices(data.shape)
    
    # Calculate the x and y coordinates of the data
    x = (X*data).sum()/total
    y = (Y*data).sum()/total          
    
    # Calculate the column and row of the data
    col = data[:, int(y)]
    row = data[int(x), :]
    
    # Calculate the sum of the column and row
    colsum = np.sum(col)
    rowsum = np.sum(row)
    
    # Calculate the width of the x and y coordinates
    width_x = np.sqrt(np.abs((np.arange(col.size)-x)**2*col).sum()/colsum)
    width_y = np.sqrt(np.abs((np.arange(row.size)-y)**2*row).sum()/rowsum)
    
    # Calculate the skew of the x and y coordinates
    skew_x = ((np.arange(col.size)-x)**3*col).sum()/colsum/(width_x**3)
    skew_y = ((np.arange(row.size)-y)**3*row).sum()/rowsum/(width_y**3)
    
    # Calculate the kurtosis of the x and y coordinates
    kurt_x = ((np.arange(col.size)-x)**4*col).sum()/colsum/(width_x**4)
    kurt_y = ((np.arange(row.size)-y)**4*row).sum()/rowsum/(width_y**4)
    
    # Return the height, x, y, width_x, width_y, skew_x, skew_y, kurt_x, kurt_y values
    return height, x, y, width_x, width_y, skew_x, skew_y, kurt_x, kurt_y    


def aper_phot(img,xc=0,yc=0):

    # flux weighted centroid
    #if xc==0 or yc==0:
    #    yc, xc = np.unravel_index(np.argmax(img,axis=None),  img.shape)
    xv,yv = mesh_box([xc,yc],4)

    # clamp bounds to image size
    xv = np.clip(xv,0,img.shape[1]-1)
    yv = np.clip(yv,0,img.shape[0]-1)

    try:
        wx = np.sum(np.unique(xv)* img[yv,xv].sum(0))/np.sum( img[yv,xv].sum(0))
        wy = np.sum(np.unique(yv)* img[yv,xv].sum(1))/np.sum( img[yv,xv].sum(1))
    except:
        print(xc,yc)
        xv,yv = mesh_box([xc,yc],1)
        wx = np.sum(np.unique(xv)* img[yv,xv].sum(0))/np.sum( img[yv,xv].sum(0))
        wy = np.sum(np.unique(yv)* img[yv,xv].sum(1))/np.sum( img[yv,xv].sum(1))
        
    # loop through aper sizes
    apers = []; bgs = []; npps = [];
    for r in np.arange(2,3,0.15):
        area, bg, npp, size = phot( img, wx, wy, r=r, dr=6)
        apers.append(area); bgs.append(bg); npps.append(npp)

    return wx, wy, apers, bgs, npps

if __name__ == '__main__':

    args = parse_args()
    args.target = ''.join(args.target.split(' '))

    dirs = glob.glob(args.datadir+"*/")

    # for each system
    for i in range(len(dirs)):
        
        sname = dirs[i].lower().split('/')[-2]

        # process specific target
        if args.target == "jwst":
            with open('jwst_targets.txt') as f:
                jwst_targets = [line.rstrip().lower() for line in f]
            if sname not in jwst_targets:
                continue
        elif args.target != 'all':
            #if args.target.lower() != dirs[i].lower():
            if args.target.lower() != os.path.basename(dirs[i][:-1]).lower():
                continue

        # alloc data used to store processed image status
        print(dirs[i])
        filelist = pickle.load(open(os.path.join(dirs[i],"filelist.pkl"),"rb"))

        if "imagetype" not in filelist:
            filelist["imagetype"] = ['none']*len(filelist["wavelength"])
        
        # alloc photometry data
        SV = {}
        waves = [filelist['wavelength'][j].strip() for j in range(len(filelist['wavelength']))]
        for f in np.unique(waves):
            
            SV[f] = {
                'aor':[],
                'time':[],
                'frame':[],
                'exptime':[],
                'flux_conv':[],
                'gain':[],
                'aper_flux':[],
                'aper_err':[],
                'aper_xcent':[],
                'aper_ycent':[],
                'aper_npp':[],
                'aper_bg':[],
                'aper_moments':[],   
            }

        # load priors (used to find star)
        if os.path.exists(dirs[i]+"centroid.json"):
            priors = json.load(open(dirs[i]+"centroid.json","r"))
        else:
            priors = {"centroid":{}}

        # reset list used to track processed images
        if args.reset:
            filelist["imagetype"] = ['none']*len(filelist["wavelength"])
            #pickle.dump(filelist,open(dirs[i]+"filelist.pkl","wb"))
        else: # load data
            if os.path.exists(dirs[i]+"photometry.pkl"):
                print('loading data...')
                SV = pickle.load(open(dirs[i]+"photometry.pkl","rb"))

        print("Images per filter:",dict(Counter(waves)))

        if args.aor:
            if args.aor not in filelist['aor']:
                print("AOR not found in filelist")
                continue
            
        # for each image
        for j in range(len(filelist['accessUrl'])):

            # select aor
            if args.aor:
                if filelist['aor'][j] != args.aor:
                    continue

            #process specific filter
            if args.filter != 'all':
                if args.filter.lower() not in filelist['wavelength'][j]:
                    continue

            # don't reprocess data unless "mode" is specified
            if filelist["imagetype"][j] != 'none': # full/sub
                # print(filelist["accessUrl"][j], filelist["imagetype"][j], "skipped")
                if args.mode == 'all':
                    continue

            # download image into memory
            try:
                url = filelist['accessUrl'][j].strip()
                response = urllib.request.urlopen(url,timeout=30)
            except:
                print(url,"failed to download", aor)
                continue

            # open fits
            try:
                hdulist = fits.open(BytesIO(response.read()))
            except:
                print(url,"failed to open")
                continue
            f = filelist['wavelength'][j].strip()

            # measure flux in each frame
            for hdu in hdulist:

                # find star in image
                aor = hdu.header.get("aorkey")
                centroid = priors.get('centroid',{}).get(str(aor),[])

                if len(centroid) == 0:
                    # query simbad for proper motion + use WCS
                    if sname in target_coords:
                        centroid = find_target(dirs[i].split('/')[1], hdu, 
                                                stats=target_coords[sname])
                    else:
                        centroid = find_target(dirs[i].split('/')[1], hdu)
                        
                        try:
                            centroid = find_target(dirs[i].split('/')[1], hdu)
                        except:
                            print(f"failed to find {sname} in {url}")
                            try:
                                centroid = find_target(hdu.header['OBJECT'], hdu)
                            except:
                                centroid = [-1,-1]

                    # if subframe mode, set centroid to center
                    if hdulist[0].header['READMODE'] == 'SUB':
                        dcube = hdulist[0].data[0]
                        # replace nans with median of image
                        dcube[np.isnan(dcube)] = np.nanmedian(dcube)
                        yc, xc = np.unravel_index(np.argmax(dcube,axis=None),  dcube.shape)

                        # make sure neighboring pixels are not background (ie cosmic ray)
                        max_mask = dcube == np.nanmax(dcube)
                        dmax_mask = binary_dilation(max_mask, iterations=2)
                        xmax_mask = np.logical_xor(max_mask, dmax_mask)
                        mean_mask = np.nanmean(dcube[xmax_mask])
                        bg_level = np.nanpercentile(dcube,75)
                        if mean_mask/bg_level < 2:                        
                            # redo we found a cosmic ray
                            # mask out cosmic ray with mean_mask
                            dcube[max_mask] = mean_mask
                            yc, xc = np.unravel_index(np.argmax(dcube,axis=None),  dcube.shape)

                        centroid = [float(xc),float(yc)]
                        print(centroid)

                    # centroid plot
                    fig,ax = plt.subplots(1)
                    hdu = hdulist[0]
                    if hdu.data.ndim == 3:
                        data = np.nanmedian(hdu.data,0)
                        ax.imshow(data,vmin=np.nanmin(data), vmax=0.25*np.nanmax(data))
                    if hdu.data.ndim == 2:
                        ax.imshow(hdu.data,vmin=np.nanmin(hdu.data), vmax=0.25*np.nanmax(hdu.data))
                    ax.plot(centroid[0], centroid[1],'rx',alpha=0.5,ms=25)
                    plt.tight_layout()
                    fig.savefig(dirs[i]+"{}_centroid.png".format(str(aor)))
                    plt.close()

                    priors['centroid'][str(aor)] = list(centroid)
                    json.dump(priors, open(dirs[i]+"centroid.json","w"))


                check = (centroid[0] > 0) & (centroid[1] > 0) & (centroid[0] < hdu.data.shape[1]) & (centroid[1] < hdu.data.shape[0]) 
                if check:
                    # todo check for star near edge of image, will fail later on if so
                    pass
                else:
                    filelist["imagetype"][j] = "bad centroid"
                    print("bad centroid",j,centroid,aor)
                    continue
                    #print(dirs[i])
                    #print('centroid:',centroid)
                    #print('data shape:',hdu.data.shape)
                    #import pdb; pdb.set_trace()
                    #hdu = hdulist[0]
                    #plt.imshow(hdu.data,vmin=np.nanmin(hdu.data.flatten()), vmax=0.25*np.nanmax(hdu.data.flatten()) )
                    #plt.plot(centroid[0], centroid[1],'rx',alpha=0.5)
                    #plt.show()

                if hdu.data.ndim == 2:
                    filelist["imagetype"][j] = "full"

                    # skip 
                    if args.mode == "sub":
                        continue
                    
                    if hdu.header.get('exptype')!='sci':
                        continue

                    start = hdu.header.get('MJD_OBS') + 2400000.5

                    data = hdu.data.copy()
                    exptime = hdu.header.get('ATIMEEND') - hdu.header.get('AINTBEG')
                    # convert from MJy/sr to e-
                    data /= float(hdu.header.get('FLUXCONV',0.1257))
                    data *= float(hdu.header.get('GAIN',3.7))
                    data *= float(exptime) # photons
                    data[np.isnan(data)] = 0
                    data[np.isinf(data)] = 0

                    try:
                        # aperture photometry
                        xc, yc = centroid
                        wx, wy, apers, bgs, npps = aper_phot(data, xc, yc)

                        # PSF photometry
                        #psf_pars = psf_phot(data, xc, yc)

                        # save data
                        SV[f]['aor'].append(aor)
                        SV[f]['time'].append(start)
                        SV[f]['frame'].append(aor)
                        SV[f]['aper_flux'].append(apers)
                        SV[f]['aper_npp'].append(npps)
                        SV[f]['aper_bg'].append(bgs)
                        SV[f]['aper_xcent'].append(wx)
                        SV[f]['aper_ycent'].append(wy)
                        #SV[f]['psf_pars'].append(psf_pars)
                        SV[f]['gain'].append(hdu.header.get('GAIN',3.7))
                        SV[f]['flux_conv'].append(hdu.header.get('FLUXCONV',0.1257))
                        SV[f]['exptime'].append(exptime)
                    except (ValueError, IndexError) as ex:
                        print(ex,)
                        # likely a cosmic ray made a bright pixel near the edge
                        # aperture routines bug out when computing npp
                        print("skipped",j)
                elif hdu.data.ndim == 3:
                    filelist["imagetype"][j] = "sub"

                    # skip 
                    if args.mode == "full":
                        continue

                    if (hdu.size != 0) and (hdu.header.get('exptype')=='sci'):
                        start = hdu.header.get('MJD_OBS') + 2400000.5

                        dcube = hdu.data.copy()
                        # convert from Mjy/sr to DN/s then to e/s and finally e
                        dcube /= float(hdu.header.get('FLUXCONV',0.1257))
                        dcube *= float(hdu.header.get('GAIN',3.7))

                        idur = hdu.header.get('ATIMEEND') - hdu.header.get('AINTBEG')
                        nimgs = dcube.shape[0]
                        exptime = idur/nimgs # sec
                        dcube[np.isnan(dcube)] = 0
                        dcube[np.isinf(dcube)] = 0
                        dcube *= float(exptime)

                        # loop through dcube
                        for k in range(dcube.shape[0]):

                            #try:
                                # aperture photometry
                                # estimate centroid as brightest pixel
                                # yc, xc = np.unravel_index(np.argmax(dcube[k],axis=None),  dcube[k].shape)

                                # # make sure neighboring pixels are not background (ie cosmic ray)
                                # max_mask = dcube[k] == np.max(dcube[k])
                                # dmax_mask = binary_dilation(max_mask, iterations=1)
                                # xmax_mask = np.logical_xor(max_mask, dmax_mask)
                                # mean_mask = dcube[k][xmax_mask].mean()
                                
                                # if mean_mask < np.percentile(dcube[k],75):
                                #     # redo we found a cosmic ray
                                #     # mask out cosmic ray with mean_mask
                                #     dcube[k][max_mask] = mean_mask
                                #     yc, xc = np.unravel_index(np.argmax(dcube[k],axis=None),  dcube[k].shape)

                                # replace nans with median of image
                                dcube[k][np.isnan(dcube[k])] = np.nanmedian(dcube[k])
                                yc, xc = np.unravel_index(np.argmax(dcube[k],axis=None),  dcube[k].shape)

                                # make sure neighboring pixels are not background (ie cosmic ray)
                                max_mask = dcube[k] == np.nanmax(dcube[k])
                                dmax_mask = binary_dilation(max_mask, iterations=2)
                                xmax_mask = np.logical_xor(max_mask, dmax_mask)
                                mean_mask = np.nanmean(dcube[k][xmax_mask])
                                bg_level = np.nanpercentile(dcube[k],75)
                                if mean_mask/bg_level < 2:                        
                                    # redo we found a cosmic ray
                                    # mask out cosmic ray with mean_mask
                                    dcube[k][max_mask] = mean_mask
                                    yc, xc = np.unravel_index(np.argmax(dcube[k],axis=None),  dcube[k].shape)

                                #centroid = [float(xc),float(yc)]

                                # if xc,yc is on edge of detector then skip
                                if xc == 0 or yc == 0 or xc == dcube[k].shape[1]-1 or yc == dcube[k].shape[0]-1:
                                    print("skipped, centroid on edge of detector",j,k)
                                    continue

                                wx, wy, apers, bgs, npps = aper_phot(dcube[k], xc, yc)

                                # extract center half of image
                                # 16x16 region around centroid
                                centerslice = dcube[k][dcube[k].shape[0]//4:-dcube[k].shape[0]//4,dcube[k].shape[1]//4:-dcube[k].shape[1]//4]
                                moments = calc_moments(centerslice)

                                # save data
                                SV[f]['aor'].append(aor)
                                SV[f]['time'].append(start+(k+1)*exptime/24/60/60)
                                SV[f]['frame'].append(k)
                                SV[f]['aper_flux'].append(apers)
                                SV[f]['aper_npp'].append(npps)
                                SV[f]['aper_bg'].append(bgs)
                                SV[f]['aper_xcent'].append(wx)
                                SV[f]['aper_ycent'].append(wy)
                                SV[f]['aper_moments'].append(moments)
                                #SV[f]['psf_pars'].append(psf_pars)
                                SV[f]['exptime'].append(exptime)
                                SV[f]['gain'].append(hdu.header.get('GAIN',3.7))
                                SV[f]['flux_conv'].append(hdu.header.get('FLUXCONV',0.1257))
                            #except (ValueError, IndexError) as e:
                                #import pdb; pdb.set_trace()
                                # likely a cosmic ray made a bright pixel near the edge
                                # aperture routines bug out when computing npp
                                #print("skipped",j,k,e)
                else:
                    pass

            del(response)
            hdulist.close()
            del(hdulist)
            
            # checkpoint save
            if j%100 == 0 and j>0:
                pickle.dump(filelist, open(dirs[i]+"filelist.pkl","wb"))
                if args.aor != 0:
                    pickle.dump(SV,open(dirs[i]+f"photometry_{args.aor}.pkl","wb"))
                else:
                    pickle.dump(SV,open(dirs[i]+"photometry.pkl","wb"))
                print(j,'saved')

        # one file save at the end if data exists
        pickle.dump(filelist, open(dirs[i]+"filelist.pkl","wb"))
        if args.aor != 0:
            pickle.dump(SV,open(dirs[i]+f"photometry_{args.aor}.pkl","wb"))
        else:
            pickle.dump(SV,open(dirs[i]+"photometry.pkl","wb"))

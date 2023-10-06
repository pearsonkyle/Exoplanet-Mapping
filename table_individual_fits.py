import os
import json
import glob
import pickle
import argparse
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

def parse_args():
    parser = argparse.ArgumentParser()

    help_ = "Choose a target to process"
    parser.add_argument("-t", "--target", help=help_, type=str, default="HD189")

    help_ = "Choose a filter (3.6, 4.5, 5.8, 8.0)"
    parser.add_argument("-f", "--filter", help=help_, type=str, default="all")

    help_ = "Choose a planet (b, c, d, ...)"
    parser.add_argument("-p", "--planet", help=help_, type=str, default="all")
    
    help_ = "Directory containing list of stars"
    parser.add_argument("-d", "--datadir", help=help_, default="DATA/", type=str)

    help_ = "Choose a frame number to process (0=all or 1-64)"
    parser.add_argument('-fr', '--frame', help=help_, default=0, type=int)

    return parser.parse_args()

def time_bin(time, flux, dt=1./(60*24)):
    bins = int(np.floor((max(time) - min(time))/dt))
    bflux = np.zeros(bins)
    btime = np.zeros(bins)
    for i in range(bins):
        mask = (time >= (min(time)+i*dt)) & (time < (min(time)+(i+1)*dt))
        if mask.sum() > 0:
            bflux[i] = np.nanmedian(flux[mask])
            btime[i] = np.nanmedian(time[mask])
    zmask = (bflux==0) | (btime==0) | np.isnan(bflux) | np.isnan(btime)
    return btime[~zmask], bflux[~zmask]

if __name__ == "__main__":
    args = parse_args()

    dirs = glob.glob(args.datadir+"*/")

    # for each system
    for i in range(len(dirs)):

        # process specific target
        sname = dirs[i].lower().split('/')[-2]
        if args.target == "jwst":
            with open('jwst_targets.txt') as f:
                jwst_targets = [line.rstrip().lower() for line in f]
            if sname not in jwst_targets:
                continue
        elif args.target != 'all':
            if args.target.lower() not in dirs[i].lower():
                continue

        print(dirs[i])
        filters = glob.glob(os.path.join(dirs[i],'IRAC*/'))

        # load data
        try:
            priors = json.load(open(dirs[i]+"prior.json","r"))
            #photometry = pickle.load(open(dirs[i]+"photometry.pkl","rb"))
        except:
            print(' no prior or no photometry files')
            continue

        # photometry{filter}.keys()
        # ['time', 'frame', 'aper_flux', 'aper_err', 'aper_xcent', 'aper_ycent', 'aper_npp', 'aper_bg', 
        # 'psf_flux', 'psf_err', 'psf_xcent', 'psf_ycent', 'psf_xsigma', 'psf_ysigma', 'psf_rot'])

        print(dirs[i])

        # alloc data
        sv = {}

        # for each filter
        for fdir in filters:
            f = os.path.basename(fdir[:-1])

            #process specific filter
            if args.filter != 'all':
                if args.filter.lower() not in f:
                    continue
 
            # collect data from the observations
            data = {
                'time':[],
                'frames':[],
                'flux':[],
                'err':[],
                'xcent':[],
                'ycent':[],
                'npp':[],
                'ramp':[],
                'detrended':[],
                'residuals':[],
                'phase':[],
                'transit':[],
                'pars':[],
                'errs':[],
            }

            no_aper_keys = ['hotspot_lon', 'hotspot_lon_std', 'hotspot_lat', 'hotspot_lat_std', 'max_amp', 'max_amp_std', 'btemp', 'btemp_std', 'photon_noise' ]
            for k in no_aper_keys:
                data[k] = []

            aors = []
            tmids = [] # used to sort results
            # skip filter if no data
            observations = glob.glob(os.path.join(fdir,f"*_eclipse_{args.frame}"))

            for obs in observations:
                if 'global' in obs:
                    continue

                # load data
                print(os.path.basename(obs))
                try:
                    photometry = pickle.load(open(os.path.join(obs,"data.pkl"),"rb"))
                except:
                    continue

                # find aor
                try:
                    aorfile = glob.glob(os.path.join(obs,"aor*txt"))
                    aor = aorfile[0].split('/')[-1].split('_')[1].split('.')[0]
                    aors.append(aor)
                except:
                    pass

                # add data to the list
                for key in data.keys():
                    if key in no_aper_keys:
                        data[key].append(photometry[key])
                    else:
                        data[key.replace("aper_","")].append(photometry["aper_"+key])
                tmids.append(photometry["aper_pars"]['tmid']+photometry["aper_pars"]['dt'])

            nobs = len(data['time'])
            tmids  = np.array(tmids)
            si = np.argsort(tmids)

            # try:
            #     # concatenate each key
            #     for key in data.keys():
            #         data[key] = np.concatenate(data[key])
            # except:
            #     # move onto next filter - usually only 1 observation
            #     continue 

            # for each observation
            for ii in range(nobs):
                i = ii #si[ii] # order data by time

                # calc phase correction 
                dphase = data['pars'][i]['dt']*2/priors['b']['period']

                # detrended
                emid = data['pars'][i]['tmid'] + 0.5*data['pars'][i]['per'] - data['pars'][i]['dt']

                bt, bd = time_bin(data['phase'][i]*priors['b']['period']+dphase, data['detrended'][i])
                bp = bt/priors['b']['period']

                # residuals
                bt, br = time_bin(data['phase'][i]*priors['b']['period'], 1e6*data['residuals'][i]/np.median(data['flux'][i]))
                stdev = np.nanstd(1e6*data['residuals'][i]/np.median(data['flux'][i]))
                stdev_bin = np.nanstd(br)

                # bin data with various samples and plot stdev vs bin size
                ppmres = 1e6*data['residuals'][i]/np.median(data['flux'][i])
                stdev_calc = [1e6*data['photon_noise'][i]]
                stdevs=[np.std(ppmres)]
                dt = np.diff(data['time'][i]).mean()
                for b in np.logspace(-2, 1, 15):
                    #bt, br = time_bin(data['phase'][i]*priors['b']['period'], ppmres, b/(24*60))
                    npts = int(b/(dt*24*60))

                    stdevs.append(ppmres[:npts*(len(ppmres)//npts)].reshape(-1, npts).mean(axis=1).std())
                    stdev_calc.append( stdev_calc[0]/np.sqrt(npts) )

                bins = np.concatenate([[dt*24*60], np.logspace(-2, 1, 15)])
                
                # plot Gaussian of hotspot
                degrees = np.linspace(-45, 45, 100)
                long_hist = np.exp( -(degrees-data['hotspot_lon'][ii])**2 / (2*data['hotspot_lon_std'][ii]**2) )
                lat_hist = np.exp( -(degrees-data['hotspot_lat'][ii])**2 / (2*data['hotspot_lat_std'][ii]**2) )
                long_hist /= np.sum(long_hist)
                lat_hist /= np.sum(lat_hist)
                
                wf = np.array(data['flux'][i]) / np.array(data['detrended'][i])

                date = Time(data['time'][i].min(),format='jd').isot.split('T')[0]

                # Table 1: AOR, DATE, X, Y, NPP, N Observations, residual stdev, photon noise, photon noise factor
                print("%s & %s & %.2f & %.2f & %.2f & %d & %d & %d & %.2f \\\\" % (aors[i], date, data['xcent'][i].mean(), data['ycent'][i].mean(), data['npp'][i].mean(), len(data['time'][i]), stdevs[0], stdev_calc[0], stdevs[0]/stdev_calc[0]))

                # spherical harmonic errors
                Y00 = "%.2f $\pm$ %.2f" % (data['pars'][i]['y00']*1e6, data['errs'][i]['y00']*1e6)
                try:
                    Y10 = "%.2f $\pm$ %.2f" % (data['pars'][i]['y10']*1e6, data['errs'][i]['y10']*1e6)
                except:
                    Y10 = "%.2f $\pm$ 0" % (data['pars'][i]['y10']*1e6)

                Y11 = "%.2f $\pm$ %.2f" % (data['pars'][i]['y11']*1e6, data['errs'][i]['y11']*1e6)
                Y1m1 = "%.2f $\pm$ %.2f" % (data['pars'][i]['y1m1']*1e6, data['errs'][i]['y1m1']*1e6)

                # calculations
                hotlong = "%.2f $\pm$ %.2f" % (data['hotspot_lon'][i], data['hotspot_lon_std'][i])
                hotlat = "%.2f $\pm$ %.2f" % (data['hotspot_lat'][i], data['hotspot_lat_std'][i])
                edepth = "%.2f $\pm$ %.2f" % (data['max_amp'][i]*1e6, data['max_amp_std'][i]*1e6)
                daytemp = "%.0f $\pm$ %.0f" % (data['btemp'][i], data['btemp_std'][i])

                # Table 2: AOR, Emid +- Err, Eclipse Depth, Hotspot Longitude, Hotspot Latitude, Daytime Temperature
                print("%s & %.5f $\pm$ %.5f & %s & %s & %s & %s \\\\" % (aors[i], emid, data['errs'][i]['dt'], edepth, hotlong, hotlat, daytemp))


                # print("Residuals Stdev [ppm]   & %.2f \\\\" % (stdevs[0]))
                # print("Photon Noise [ppm]  & %.2f \\\\" % (stdev_calc[0]))
                # print("Photon Noise Factor  & %.2f \\\\" % (stdevs[0]/stdev_calc[0]))
                # print("N. Observations  & %i \\\\" % (len(data['time'][i])))
                # print("Average X-Centroid [pix]  & %.2f $\pm$ %.2f \\\\" % (np.nanmean(data['xcent'][i]), np.nanstd(data['xcent'][i])))
                # print("Average Y-Centroid [pix]  & %.2f $\pm$ %.2f \\\\" % (np.nanmean(data['ycent'][i]), np.nanstd(data['ycent'][i])))
                # print("Average Noise Pixel  & %.2f $\pm$ %.2f \\\\" % (np.nanmean(data['npp'][i]), np.nanstd(data['npp'][i])))

import os
import json
import glob
import pickle
import argparse
import numpy as np
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

            no_aper_keys = ['hotspot_lon', 'hotspot_lon_std', 'hotspot_lat', 'hotspot_lat_std', 'max_amp', 'max_amp_std', 'btemp', 'btemp_std' ]
            for k in no_aper_keys:
                data[k] = []

            aors = []
            tmids = [] # used to sort results
            # skip filter if no data
            observations = glob.glob(os.path.join(fdir,f"*_eclipse_{args.frame}"))

            for obs in observations:
                if 'global' not in obs:
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

            fig = plt.figure(figsize=(7,5))
            ax = plt.subplot2grid((3,1),(0,0),rowspan=2)
            axr = plt.subplot2grid((3,1),(2,0))
            plt.subplots_adjust(hspace=0.035)

            title = f"{sname.upper()} b {f}"
            si = np.argsort(data['phase'][0])
            # for each observation
            for i in range(len(data['time'])):
                #ax.plot(data['phase'][i], data['detrended'][i], 'k.',ms=1, alpha=0.1)
                ax.set_xlim([0.47,0.53])
                axr.set_xlim([0.47,0.53])
                mask = (data['phase'][i][si] > 0.47) & (data['phase'][i][si] < 0.53)
                bt, bd = time_bin(data['phase'][i][si][mask]*data['pars'][0]['per'], data['detrended'][i][si][mask], dt=1./(60*24))
                bp = bt/data['pars'][0]['per']
                ax.set_ylabel('Rel. Flux', fontsize=14)
                ax.set_title(title, fontsize=14)
                ax.set_ylim([0.9995,1.0025])
                ax.set_xticklabels([])
                ax.grid(True,ls='--',alpha=0.5)
                
                # residuals
                stdev = np.nanstd(1e6*data['residuals'][i][mask]/np.median(data['flux'][i]))
                bt, br = time_bin(data['phase'][i][si][mask]*priors['b']['period'], data['residuals'][i][si][mask]/np.median(data['flux'][i]))
                stdev_bin = np.nanstd(br)

                # sigma clip residuals
                smask = np.abs(br) < 3*stdev_bin

                # eclipse plot
                ax.plot(data['phase'][i][si], data['transit'][i][si], 'r-', lw=2, zorder=2,alpha=0.75)
                ax.errorbar(bp[smask], bd[smask], yerr=stdev_bin, ls='none',marker='.',color='black', alpha=0.5, zorder=1)

                # smooth binned data with savgol
                sbd = savgol_filter(bd[smask], 11, 3)


                # residuals plot
                bp = bt/data['pars'][0]['per']
                axr.errorbar(bp[smask], br[smask]*1e6, yerr=stdev_bin*1e6, label=f"stdev: {stdev_bin*1e6:.1f} ppm", ls='none',marker='.',color='black', alpha=0.5, zorder=1)
                axr.set_ylabel("Res. [ppm]", fontsize=14)
                axr.legend(fontsize=12, loc='upper left')
                
                axr.grid(True,ls='--',alpha=0.5)
                axr.set_xlabel('Phase', fontsize=14)


                # bin data with various samples and plot stdev vs bin size
                ppmres = 1e6*data['residuals'][i]/np.median(data['flux'][i])
                stdevs=[np.std(ppmres)]
                dt = np.diff(data['time'][i]).mean()
                for b in np.logspace(-2, 1, 15):
                    #bt, br = time_bin(data['phase'][i]*priors['b']['period'], ppmres, b/(24*60))
                    npts = int(b/(dt*24*60))

                    stdevs.append(ppmres[:npts*(len(ppmres)//npts)].reshape(-1, npts).mean(axis=1).std())

                bins = np.concatenate([[dt*24*60], np.logspace(-2, 1, 15)])

                # plot Gaussian of hotspot
                degrees = np.linspace(-45, 45, 100)
                long_hist = np.exp( -(degrees-data['hotspot_lon'][i])**2 / (2*data['hotspot_lon_std'][i]**2) )
                lat_hist = np.exp( -(degrees-data['hotspot_lat'][i])**2 / (2*data['hotspot_lat_std'][i]**2) )
                long_hist /= np.sum(long_hist)
                lat_hist /= np.sum(lat_hist)

                wf = np.array(data['flux'][i]) / np.array(data['detrended'][i])

                emid = data['pars'][i]['tmid'] + 0.5*data['pars'][i]['per'] - data['pars'][i]['dt']

                try:
                    print("E$_{mid}$ [BJD]  & Mid Eclipse Time   & %.5f $\pm$ %.5f \\\\" % (emid, data['errs'][i]['dt']))
                except:
                    pass

                try:
                    print("AOR & %s \\\\" % aors[i])
                except:
                    pass

                # spherical harmonic errors
                print("Y$_00$ [ppm]   & Uniform Brightness   & %.2f $\pm$ %.2f \\\\" % (data['pars'][i]['y00']*1e6, data['errs'][i]['y00']*1e6))
                try:
                    print("Y$_10$ [ppm]   & Day-Night Amplitude   & %.2f $\pm$ %.2f \\\\" % (data['pars'][i]['y10']*1e6, data['errs'][i]['y10']*1e6))
                except:
                    print("Y$_10$ [ppm]   & Day-Night Amplitude   & %.2f $\pm$ 0 \\\\" % (data['pars'][i]['y10']*1e6))

                try:
                    print("Y$_11$ [ppm]   & Hot-spot Offset   & %.2f $\pm$ %.2f \\\\" % (data['pars'][i]['y11']*1e6, data['errs'][i]['y11']*1e6))
                except:
                    pass
                try:
                    print("Y$_1m1$ [ppm]   & Hot-spot Latitude  & %.2f $\pm$ %.2f \\\\" % (data['pars'][i]['y1m1']*1e6, data['errs'][i]['y1m1']*1e6))
                except:
                    pass

                # hot-spot longitude
                try:
                    print("Hotspot Longitude [deg]   & %.2f $\pm$ %.2f \\\\" % (data['hotspot_lon'][i], data['hotspot_lon_std'][i]))
                    print("Hotspot Latitude [deg]   & %.2f $\pm$ %.2f \\\\" % (data['hotspot_lat'][i], data['hotspot_lat_std'][i]))
                    print("Eclipse Depth [ppm]   & %.2f $\pm$ %.2f \\\\" % (data['max_amp'][i]*1e6, data['max_amp_std'][i]*1e6))
                    print("Day-side Temperature [K]   & %.0f $\pm$ %.0f \\\\" % (data['btemp'][i], data['btemp_std'][i]))
                except:
                    pass

                print("Residuals Stdev [ppm]   & %.2f \\\\" % (stdevs[0]))
                print("Binned Residual Stdev [ppm]   & %.2f \\\\" % (stdev_bin*1e6))
                print("N. Observations  & %i \\\\" % (len(data['time'][i])))

            plt.tight_layout()
            plt.savefig(f"{title.replace('.','').replace(' ','_')}_global.png", dpi=300)
            plt.show()

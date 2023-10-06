import os
import json
import glob
import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
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
                aorfile = glob.glob(os.path.join(obs,"aor*txt"))
                aor = aorfile[0].split('/')[-1].split('_')[1].split('.')[0]
                aors.append(aor)

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

            # loop through planets

            fig, ax = plt.subplots(nobs, 7, figsize=(14, 1.7*nobs))
            fig.subplots_adjust(hspace=0.08, wspace=0.25, left=0.06, right=0.98, top=0.94, bottom=0.034)
            fig.suptitle(f"{sname.upper()} b IRAC {args.filter}" +r" $\mu$m Eclipse", fontdict=dict(weight='bold'))

            mins = []
            maxs = []
            for flux in data['flux']:
                mins.append(np.nanmin(flux))
                maxs.append(np.nanmax(flux))
            minf = np.nanmin(mins)
            maxf = np.nanmax(maxs)

            # for each observation
            for ii in range(nobs):
                print(ii)
                i = si[ii] # order data by time
                ax[ii,0].plot(data['phase'][i], data['flux'][i], 'k.',ms=1, alpha=0.25)
                #ax[ii,0].set_ylim([minf, maxf])
                ax[ii,0].set_xlim([0.45,0.55])

                # xcent,ycent
                ax[ii,2].plot(data['phase'][i], data['xcent'][i], 'r^', alpha=0.25, ms=1)#, label="X - red")
                ax[ii, 2].text(0.55, np.nanmean(data['xcent'][i]), ' X', ha='left',
                               va='center', c= 'r')
                ax[ii,2].plot(data['phase'][i], data['ycent'][i], 'b.', alpha=0.25, ms=1)#, label="Y - blue")
                ax[ii, 2].text(0.55, np.nanmean(data['ycent'][i]), ' Y', ha='left',
                               va='center', c='b')
                ax[ii,2].grid(True, ls='--')
                ax[ii,2].set_xlim([0.45,0.55])

                # npp
                ax[ii,1].plot(data['phase'][i], data['npp'][i], 'k.',ms=1, alpha=0.25)
                ax[ii,1].grid(True, ls='--')
                ax[ii,1].set_xlim([0.45,0.55])

                # calc phase correction 
                dphase = data['pars'][i]['dt']*2/priors['b']['period']

                # detrended
                emid = data['pars'][i]['tmid'] + 0.5*data['pars'][i]['per'] - data['pars'][i]['dt']

                bt, bd = time_bin(data['phase'][i]*priors['b']['period']+dphase, data['detrended'][i])
                bp = bt/priors['b']['period']
                ax[ii,4].plot(bp+dphase*0.5, bd, 'k.', alpha=0.5)#, label=f"AOR {aors[i]}")
                ax[ii,4].plot(data['phase'][i]+dphase, data['transit'][i], 'r-')
                ax[ii,4].set_ylim([0.9995, 1.003])
                ax[ii,4].grid(True, ls='--')
                ax[ii,4].set_xlim([0.45,0.55])

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
                # ax[ii,5].semilogx(dt*24*60, stdev_calc[0], 'r-',lw=3)#, label=r" $\sigma_{photon}$ ="+"\n"+f"  {stdev_calc[0]:.0f} ppm")

                ax[ii,5].text(10 ** -3, 10**2, r" $\sigma_{photon}$ ="+"\n"+f"  {stdev_calc[0]:.0f} ppm", ha='left', va='top', c='r')


                ax[ii,5].semilogx(dt*24*60, stdevs[0], 'k-', lw=0.5)#, label=r" $\sigma_{residuals}$ ="+"\n"+f"  {stdevs[0]:.0f} ppm")

                ax[ii, 5].text(10 ** 1, 10**3.5, r"$\sigma_{residuals}$ ="+"\n"+f"{stdevs[0]:.0f} ppm ", ha='right', va='top', alpha=0.5)

                ax[ii,5].loglog(bins, stdevs, 'k-.', alpha=0.5)
                ax[ii,5].loglog(bins, stdev_calc, 'r-', alpha=1)
                # ax[ii,5].axvline(1, ls='--', c='k')
                # ax[ii,5].legend(loc = 'lower left', fontsize=8)
                # ax[ii,5].set_xlim()
                ax[ii,5].set_xlim([10 ** -3., 10 ** 1])
                ax[ii,5].set_ylim([10**1.4,10**3.5])

                # plot Gaussian of hotspot
                degrees = np.linspace(-45, 45, 100)
                long_hist = np.exp( -(degrees-data['hotspot_lon'][ii])**2 / (2*data['hotspot_lon_std'][ii]**2) )
                lat_hist = np.exp( -(degrees-data['hotspot_lat'][ii])**2 / (2*data['hotspot_lat_std'][ii]**2) )
                # normalize
                long_hist /= np.sum(long_hist)
                lat_hist /= np.sum(lat_hist)
                ax[ii,6].plot(degrees, lat_hist, 'b-')#, label="Lat =\n"+f"  {data['hotspot_lat'][ii]:.1f}"+r"$^{\circ}\pm$"+ f"{data['hotspot_lat_std'][ii]:.1f}"+r"$^{\circ}$")
                ax[ii, 6].plot(degrees, long_hist, 'r-.')#, label=f"Lon {data['hotspot_lon'][ii]:.1f}" + r"$^{\circ}\pm$" + f"{data['hotspot_lon_std'][ii]:.1f}" + r"$^{\circ}$")

                txt = ax[ii, 6].text(data['hotspot_lat'][ii], np.nanmax(lat_hist),
                                f"Lat: {data['hotspot_lat'][ii]:.1f}"+r"$^{\circ}\pm$"+ f"{data['hotspot_lat_std'][ii]:.1f}"+r"$^{\circ}$",
                                ha='center',va='bottom',fontsize=8,c='b',weight='bold')

                txt.set_path_effects([PathEffects.withStroke(linewidth=1.5, foreground='w')])

                txt = ax[ii, 6].text(data['hotspot_lon'][ii], np.nanmax(long_hist),
                               f"Lon: {data['hotspot_lon'][ii]:.1f}" + r"$^{\circ}\pm$" + f"{data['hotspot_lon_std'][ii]:.1f}" + r"$^{\circ}$",
                               ha='center', va='bottom', fontsize=8,c='r',weight='bold')

                txt.set_path_effects([PathEffects.withStroke(linewidth=1.5, foreground='w')])

                # remove top and right spines
                ax[ii,6].spines['top'].set_visible(False)
                ax[ii,6].spines['right'].set_visible(False)
                ax[ii,6].spines['left'].set_visible(False)


                # ax[ii,6].legend(loc = 'best', fontsize=8)
                # ax[ii, 6].legend(loc='upper right', fontsize=8)
                ax[ii,6].set_yticks([])
                #ax[ii,5].set_ylim([0, 0.5])
                ax[ii, 6].set_xlim([-45, 45])

                wf = np.array(data['flux'][i]) / np.array(data['detrended'][i])
                ax[ii,3].scatter(data['xcent'][i], data['ycent'][i],c=wf/np.median(wf),cmap='jet',s=data['npp'][i]/2)
                ax[ii,3].axis('off')

                ax[ii, 0].set_ylabel(f"AOR {aors[i]}", fontdict=dict(weight='bold'))

                if ii == 0:
                    # ax[ii,2].legend(loc = 'best')
                    ax[ii,0].set_title("Flux [MJy/sr]")
                    ax[ii,2].set_title("Centroid [pix]")
                    ax[ii,1].set_title("Noise Pixel")
                    ax[ii,4].set_title("Rel. Flux")
                    ax[ii,5].set_title("Residual Stdev [ppm]")
                    ax[ii,6].set_title(r"Hotspot Location [$^\circ$]")
                    ax[ii,3].set_title("Pixel Map [X,Y]")
                    ax[-1,0].set_xlabel("Phase")
                    ax[-1,1].set_xlabel("Phase")
                    ax[-1,2].set_xlabel("Phase")
                    ax[-1,3].set_xlabel("Phase")
                    ax[-1,4].set_xlabel("Phase")
                    ax[-1,5].set_xlabel("Bin Size [min]")
                    ax[-1,6].set_xlabel(r"Offset [$^\circ$]")

                try:
                    print("E$_{mid}$ [BJD]  & Mid Eclipse Time   & %.5f $\pm$ %.5f \\\\" % (emid, data['errs'][i]['dt']))
                except:
                    pass

                print("AOR & %s \\\\" % aors[i])

                # spherical harmonic errors
                print("Y$_00$ [ppm]   & Uniform Brightness   & %.2f $\pm$ %.2f \\\\" % (data['pars'][i]['y00']*1e6, data['errs'][i]['y00']*1e6))
                try:
                    print("Y$_10$ [ppm]   & Day-Night Amplitude   & %.2f $\pm$ %.2f \\\\" % (data['pars'][i]['y10']*1e6, data['errs'][i]['y10']*1e6))
                except:
                    pass

                try:
                    print("Y$_11$ [ppm]   & Hot-spot Offset   & %.2f $\pm$ %.2f \\\\" % (data['pars'][i]['y11']*1e6, data['errs'][i]['y11']*1e6))
                except:
                    pass
                try:
                    print("Y$_1m1$ [ppm]   & Hot-spot Latitude  & %.2f $\pm$ %.2f \\\\" % (data['pars'][i]['y1m1']*1e6, data['errs'][i]['y1m1']*1e6))
                except:
                    pass

                # phase curve stuff
                # hot-spot longitude
                try:
                    print("Hotspot Longitude [deg]   & %.2f $\pm$ %.2f \\\\" % (data['hotspot_lon'][i], data['hotspot_lon_std'][i]))
                    print("Hotspot Latitude [deg]   & %.2f $\pm$ %.2f \\\\" % (data['hotspot_lat'][i], data['hotspot_lat_std'][i]))
                    print("Eclipse Depth [ppm]   & %.2f $\pm$ %.2f \\\\" % (data['max_amp'][i]*1e6, data['max_amp_std'][i]*1e6))
                    print("Day-side Temperature [K]   & %.0f $\pm$ %.0f \\\\" % (data['btemp'][i], data['btemp_std'][i]))
                except:
                    pass

                print("Residuals Stdev [ppm]   & %.2f \\\\" % (stdevs[0]))
                print("Photon Noise [ppm]  & %.2f \\\\" % (stdev_calc[0]))
                print("Photon Noise Factor  & %.2f \\\\" % (stdevs[0]/stdev_calc[0]))
                print("N. Observations  & %i \\\\" % (len(data['time'][i])))
                print("Average X-Centroid [pix]  & %.2f $\pm$ %.2f \\\\" % (np.nanmean(data['xcent'][i]), np.nanstd(data['xcent'][i])))
                print("Average Y-Centroid [pix]  & %.2f $\pm$ %.2f \\\\" % (np.nanmean(data['ycent'][i]), np.nanstd(data['ycent'][i])))
                print("Average Noise Pixel  & %.2f $\pm$ %.2f \\\\" % (np.nanmean(data['npp'][i]), np.nanstd(data['npp'][i])))

                # phase curve stuff
                # hot-spot amplitude
                #print("Hotspot Intensity [ppm]   & %.2f $\pm$ %.2f \\\\" % (sv['hotspot_amp']*1e6, sv['hotspot_amp_std']*1e6))
                # cold-spot amplitude
                #print("Coldspot Intensity [ppm]   & %.2f $\pm$ %.2f \\\\" % (sv['coldspot_amp']*1e6, sv['coldspot_amp_std']*1e6))

                # resstd = sv['aper_residuals'].std() / sv['aper_flux'].mean()
                # print("&Residuals [ppm]   & %.2f  \\\\" % (resstd*1e6))
                # print("&Photon noise [ppm]   & %.2f \\\\" % (tpars['photon_noise']*1e6))
                # print("&Noise Factor  & %.2f \\\\" % (resstd/tpars['photon_noise']))

            plt.tight_layout()

            print("\n\n*** Saving plot. ***")
            # plt.savefig("FIGURES/"+args.target + "_irac" + args.filter + ".pdf", bbox_inches='tight') # the PDF is too giant...let's comment it out for now
            plt.savefig("FIGURES/"+args.target + "_irac" + args.filter + ".png", bbox_inches='tight', dpi=300)

            plt.close()

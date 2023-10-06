import os
import json
import glob
import pickle
import argparse
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

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


# average data into bins of dt from start to finish
def time_bin(time, flux, dt=1. / (60 * 24)):
    bins = int(np.floor((max(time) - min(time)) / dt))
    bflux = np.zeros(bins)
    btime = np.zeros(bins)
    bstds = np.zeros(bins)
    for i in tqdm(range(bins)):
        mask = (time >= (min(time) + i * dt)) & (time < (min(time) + (i + 1) * dt))
        if mask.sum() > 0:
            bflux[i] = np.nanmean(flux[mask])
            btime[i] = np.nanmean(time[mask])
            bstds[i] = np.nanstd(flux[mask]) / (mask.sum() ** 0.5)
    zmask = (bflux == 0) | (btime == 0) #| np.isnan(bflux) | np.isnan(btime)
    # set all zero values to nan
    bflux[zmask] = np.nan
    return btime, bflux, bstds


if __name__ == "__main__":
    args = parse_args()

    # find directories with data in them
    dirs = glob.glob(args.datadir+"*/")

    # for each system
    for i in range(len(dirs)):

        # skip directories that don't match target
        sname = dirs[i].lower().split('/')[-2]
        if args.target.lower() not in dirs[i].lower():
            continue

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

        # find filters
        filters = glob.glob(os.path.join(dirs[i],'IRAC*/'))[::-1]

        if len(filters) <= 1:
            continue

        # alloc data
        sv = {}

        fig, ax = plt.subplots(2,2, figsize=(12,8))

        # for each filter extract data
        for fi, fdir in enumerate(filters):
            f = os.path.basename(fdir[:-1])

            #process specific filter
            if args.filter != 'all':
                if args.filter.lower() not in f:
                    continue

            # skip filter if no data
            observations = glob.glob(os.path.join(fdir,f"*_phasecurve_{args.frame}"))
            for obs in observations:
                if 'global' in obs:
                    continue

                # load data
                print(os.path.basename(obs))
                photometry = pickle.load(open(os.path.join(obs,"data.pkl"),"rb"))

                # dict_keys(['hotspot_lon', 'hotspot_lat', 'hotspot_lon_std', 'hotspot_lat_std', 'hotspot_amp', 
                # 'hotspot_amp_std', 'coldspot_amp', 'coldspot_amp_std', 'max_amp', 'max_amp_std', 'btemp', 
                # 'btemp_std', 'ntemp', 'ntemp_std', 'aper_time', 'aper_flux', 'aper_err', 'aper_xcent', 
                # 'aper_ycent', 'aper_npp', 'aper_results', 'aper_wf', 'aper_model', 'aper_transit', 
                # 'aper_residuals', 'aper_detrended', 'aper_phase', 'aper_pars', 'aper_errs', 
                # 'aper_frames', 'exp_time', 'aor', 'lat_hist', 'lon_hist', 'aper_ramp', 
                # 'noise_factor', 'photon_noise', 'mutual_info_regression'])

                mod_phase = photometry['aper_phase']%1
                si = np.argsort(mod_phase)
                # plot detrended data
                ax[0, fi].scatter(photometry['aper_phase']%1, photometry['aper_detrended'], marker='.', label=f, alpha=0.5, zorder=1, s=2, color='k')
                ax[0, fi].plot(mod_phase[si], photometry['aper_transit'][si], label="Best fit", alpha=0.5, zorder=2, lw=2, color='r')
                
                # time bin phase sorted data
                si = np.argsort(photometry['aper_phase'])
                btime, bflux, bstds = time_bin(photometry['aper_phase'][si], photometry['aper_detrended'][si], dt=10./60/24)
                ax[1, fi].errorbar(btime%1, bflux, yerr=bstds, fmt='.', label=f, alpha=0.5, zorder=1, ls='none', color='b')
                ax[1, fi].plot(mod_phase[si], photometry['aper_transit'][si], label="Best fit", alpha=0.5, zorder=2, lw=2, color='r')

                # set y-lim to be the same for all plots
                ax[0, fi].set_ylim([0.97,1.02])
                ax[1, fi].set_ylim([0.999,1.006])
                ax[0, fi].set_xlim([0,1])
                ax[1, fi].set_xlim([0,1])
                ax[0, fi].legend(loc='best')
                ax[1, fi].legend(loc='best')
                ax[0, fi].set_xlabel('Phase')
                ax[1, fi].set_xlabel('Phase')
                ax[0, fi].set_ylabel('Relative Flux')
                ax[1, fi].set_ylabel('Relative Flux')
                ax[0, fi].grid(True, ls='--')
                ax[1, fi].grid(True, ls='--')

        fig.suptitle(f"{args.target} b Phase Curve Comparison")
        plt.tight_layout()
        plt.savefig(os.path.join(dirs[i],f"phasecurve_comparison_{args.frame}.png"))
        plt.close()


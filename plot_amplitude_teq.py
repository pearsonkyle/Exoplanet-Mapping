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
    parser.add_argument("-t", "--target", help=help_, type=str, default="")

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

    # create figure
    fig, ax = plt.subplots(1, figsize=(12,8))

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
                phasecurve = lambda x: x # accidently pickled a function called phasecurve
                photometry = pickle.load(open(os.path.join(obs,"data.pkl"),"rb"))

                # dict_keys(['hotspot_lon', 'hotspot_lat', 'hotspot_lon_std', 'hotspot_lat_std', 'hotspot_amp', 
                # 'hotspot_amp_std', 'coldspot_amp', 'coldspot_amp_std', 'max_amp', 'max_amp_std', 'btemp', 
                # 'btemp_std', 'ntemp', 'ntemp_std', 'aper_time', 'aper_flux', 'aper_err', 'aper_xcent', 
                # 'aper_ycent', 'aper_npp', 'aper_results', 'aper_wf', 'aper_model', 'aper_transit', 
                # 'aper_residuals', 'aper_detrended', 'aper_phase', 'aper_pars', 'aper_errs', 
                # 'aper_frames', 'exp_time', 'aor', 'lat_hist', 'lon_hist', 'aper_ramp', 
                # 'noise_factor', 'photon_noise', 'mutual_info_regression'])
                ars = priors['b']['ars']
                Teff = priors['T*']
                Te = Teff * (1 - 0.1)**0.25 * np.sqrt(0.5/ars)
                print(f"Teq = {Te:.2f} K")

                ax.scatter(photometry['btemp'], photometry['hotspot_amp'], marker='o', color='k', alpha=0.5, s=10)

    ax.set_xlabel(r'$T_{eq}$ [K]')
    ax.set_ylabel(r'Hotspot Amplitude [ppm]')
    plt.tight_layout()
    plt.savefig(f"amplitude_teq_comparison.png")
    plt.close()


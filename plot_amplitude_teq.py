import os
import json
import glob
import pickle
import argparse
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

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

    # alloc data frame for table
    df = {'target_name':[], 'filter':[],
          'hotspot_amp':[], 'hotspot_amp_std':[], 
          'hotspot_lon':[], 'hotspot_lon_std':[], 
          'hotspot_lon_phase':[], 'hotspot_lon_phase_std':[],
          'hotspot_lat':[], 'hotspot_lat_std':[],
          'day_amp':[], 'day_amp_std':[],
          'night_amp':[], 'night_amp_std':[],
          'noise_factor':[],
    }

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
            #sv = pickle.load(open(dirs[i]+"photometry.pkl","rb"))
        except:
            print(' no prior or no sv files')
            continue

        # sv{filter}.keys()
        # ['time', 'frame', 'aper_flux', 'aper_err', 'aper_xcent', 'aper_ycent', 'aper_npp', 'aper_bg', 
        # 'psf_flux', 'psf_err', 'psf_xcent', 'psf_ycent', 'psf_xsigma', 'psf_ysigma', 'psf_rot'])

        print(dirs[i])

        # find filters
        filters = glob.glob(os.path.join(dirs[i],'IRAC*/'))[::-1]

        if len(filters) <= 1:
            continue

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
                sv = pickle.load(open(os.path.join(obs,"data.pkl"),"rb"))

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

                # save data to table
                if sv['hotspot_lon_std'] < 0.1:
                    continue

                # percent of data below 1
                mask = (sv['aper_phase'] > 0.9) & (sv['aper_phase'] < 1.1) & (sv['aper_transit'] < 1)

                # estimate nightside amplitude
                min_mask = (sv['aper_phase'] < sv['aper_phase'][mask].min()-0.01) & (sv['aper_phase'] > sv['aper_phase'][mask].min()-0.02)
                max_mask = (sv['aper_phase'] > sv['aper_phase'][mask].max()+0.01) & (sv['aper_phase'] < sv['aper_phase'][mask].max()+0.02)
                # combine masks
                combined_mask = min_mask | max_mask
                night_amp = np.mean(sv['aper_transit'][combined_mask])

                # compute dayside amplitude between phase of 0.4 and 0.6, 1.4 and 1.6
                mask = ((sv['aper_phase'] > 0.4) & (sv['aper_phase'] < 0.6) | (sv['aper_phase'] > 1.4) & (sv['aper_phase'] < 1.6)) & (sv['aper_transit'] == 1)
                min_mask = (sv['aper_phase'] < sv['aper_phase'][mask].min()-0.01) & (sv['aper_phase'] > sv['aper_phase'][mask].min()-0.02)
                max_mask = (sv['aper_phase'] > sv['aper_phase'][mask].max()+0.01) & (sv['aper_phase'] < sv['aper_phase'][mask].max()+0.02)
                # combine masks
                combined_mask = min_mask | max_mask
                day_amp = np.median(sv['aper_transit'][combined_mask])

                ax.errorbar(Te, (day_amp-night_amp)*1e6, yerr=sv['hotspot_amp_std']*1e6, fmt='o', label=f'{sname} {f} {sv["hotspot_lon"]:.1f} deg')

                # fig, ax = plt.subplots(1, figsize=(12,8))
                # ax.plot(sv['aper_phase'], sv['aper_transit'], color='black', lw=2, label='phasecurve')
                # ax.axhline(night_amp, color='cyan', ls='--', label=f'night_amp = {(night_amp-1)*1e6:.0f} ppm')
                # ax.axhline(0.5*(day_amp+1+sv['max_amp']), color='orange', ls='--', label=f'day_amp = {(day_amp-1)*1e6:.0f} ppm')
                # ax.axvline(1 - (sv['hotspot_lon'] + 180) / 360, color='r', ls='--')
                # ax.axvline(1 - (sv['hotspot_lon'] + 180) / 360 + 1, color='r', ls='--', label=f'hotspot = {sv["max_amp"]*1e6:.0f} ppm, {sv["hotspot_lon"]:.1f} deg')
                # ax.grid(True,ls='--')
                # ax.legend()
                # ax.set_xlabel('Phase')
                # ax.set_ylabel('Rel. Flux')
                # ax.set_xlim([min(sv['aper_phase']), max(sv['aper_phase'])])
                # plt.show()

                # save data to table
                df['target_name'].append(sname)
                df['filter'].append(f)
                df['hotspot_amp'].append(sv['hotspot_amp'])
                df['hotspot_amp_std'].append(sv['hotspot_amp_std'])
                df['hotspot_lon'].append(sv['hotspot_lon'])
                df['hotspot_lon_std'].append(sv['hotspot_lon_std'])
                df['hotspot_lon_phase'].append(1 - (sv['hotspot_lon'] + 180) / 360)
                df['hotspot_lon_phase_std'].append(sv['hotspot_lon_std'] / 360)
                df['hotspot_lat'].append(sv['hotspot_lat'])
                df['hotspot_lat_std'].append(sv['hotspot_lat_std'])
                df['day_amp'].append(day_amp-1)
                df['day_amp_std'].append(np.sqrt(sv['aper_errs']['y00']**2 + sv['aper_errs']['y10']**2))
                df['night_amp'].append(night_amp-1)
                df['night_amp_std'].append(np.sqrt(sv['aper_errs']['y00']**2 + sv['aper_errs']['y10']**2 + sv['aper_errs']['y11']**2 + sv['aper_errs']['y1m1']**2))
                df['noise_factor'].append(sv['noise_factor'])

    # convert to dataframe and save as csv
    import pandas as pd

    df = pd.DataFrame(df, index=None)
    df.to_csv("spitzer_phasecurves.csv", index=False)

    ax.set_xlabel(r'$T_{eq}$ [K]')
    ax.set_ylabel(r'Hotspot Amplitude [ppm]')
    plt.tight_layout()
    plt.savefig(f"amplitude_teq_comparison.png")
    print(f"Saved to amplitude_teq_comparison.png")
    plt.close()


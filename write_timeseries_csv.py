import os
import json
import glob
import pickle
import argparse
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd
from exotic.api.elca import transit
from pylightcurve import exotethys

def parse_args():
    parser = argparse.ArgumentParser()

    help_ = "Choose a target to process"
    parser.add_argument("-t", "--target", help=help_, type=str, default="WASP-52")

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

        print(dirs[i])

        # find filters
        filters = glob.glob(os.path.join(dirs[i],'IRAC*/'))[::-1]

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
                photometry = pickle.load(open(os.path.join(obs,"data.pkl"),"rb"))

                priorf = json.load(open(os.path.join(dirs[i], 'prior.json'), 'r'))
                prior = {
                    'rprs':0.1169,
                    'ars': priorf['b']['ars'],
                    'per': priorf['b']['period'],
                    'tmid': priorf['b']['t0'],
                    'inc': priorf['b']['inc'],
                    'ecc': priorf['b']['ecc'],
                    'omega': priorf['b']['omega'],
                    'u1':0.1,
                    'u2':0.1,
                    'u3':0.1,
                    'u4':0.1,
                    'T*': priorf['T*'],
                    'FE/H': priorf['FEH*'],
                    'LOGG': priorf['LOGG*']
                }
                get_prior = lambda key: float(prior[key])
                u0,u1,u2,u3 = exotethys(get_prior('LOGG'), get_prior('T*'), get_prior('FE/H'), 'TESS', method='claret', stellar_model='phoenix')
                prior['u0'],prior['u1'],prior['u2'],prior['u3'] = u0,u1,u2,u3

                fig,ax = plt.subplots(1, figsize=(10,5))
                ax.plot(photometry['aper_time'], photometry['aper_detrended'], 'k.', label='Detrended Flux')
                ax.plot(photometry['aper_time'], photometry['aper_transit'], 'r-', label='Model')

                # create a transit model
                prior['tmid'] = photometry['aper_pars']['tmid'] - photometry['aper_pars']['dtt']
                print(f"tmid: {prior['tmid']+prior['per']} +- {photometry['aper_errs']['dtt']}")
                tmodel = transit(photometry['aper_time'], prior)
                ax.plot(photometry['aper_time'], tmodel, 'g-', label='Transit Model with prior tmid')
                ax.axvline(prior['tmid']+prior['per'], color='g', linestyle='--')

                prior['tmid'] = photometry['aper_pars']['tmid'] +10*photometry['aper_errs']['dtt'] + photometry['aper_pars']['per']*0.5
                print(f"emid: {prior['tmid']} +- {photometry['aper_errs']['dt']}")
                prior['rprs'] = photometry['aper_pars']['y00']**0.5
                emodel = transit(photometry['aper_time'], prior)
                ax.plot(photometry['aper_time'], emodel, 'c-', label='Transit Model with prior emid')
                ax.axvline(prior['tmid'], color='c', linestyle='--')

                ax.legend(loc='best')
                plt.show()

                # write the timeseries to file
                df = pd.DataFrame({
                    'time': photometry['aper_time'],
                    'flux': photometry['aper_detrended'],
                    'err': photometry['aper_err']/photometry['aper_flux'],
                    'xcent': photometry['aper_xcent'],
                    'ycent': photometry['aper_ycent'],
                    'npp': photometry['aper_npp'],
                    'phase': photometry['aper_phase'],
                    'raw_flux': photometry['aper_flux'],
                    'phasecurve': photometry['aper_transit']
                })

                # Out[9]: dict_keys(['hotspot_lon', 'hotspot_lat', 'hotspot_lon_std', 
                # 'hotspot_lat_std', 'hotspot_amp', 'hotspot_amp_std', 'coldspot_amp', 
                # 'coldspot_amp_std', 'max_amp', 'max_amp_std', 'btemp', 'btemp_std', 
                # 'ntemp', 'ntemp_std', 'aper_time', 'aper_flux', 'aper_err', 'aper_xcent', 
                # 'aper_ycent', 'aper_npp', 'aper_results', 'aper_wf', 'aper_model', 
                # 'aper_transit', 'aper_residuals', 'aper_detrended', 'aper_phase', 
                # 'aper_pars', 'aper_errs', 'aper_frames', 'exp_time', 'aor', 'lat_hist', 
                # 'lon_hist', 'aper_priors', 'aper_ramp', 'aper_emid', 'noise_factor', 
                # 'photon_noise', 'mutual_info_regression'])

                # add tmid and emid to header with uncertainties
                file_header = f"# tmid: {photometry['aper_pars']['tmid']:.6f} +- {photometry['aper_errs']['dtt']:.6f}\n"
                file_header += f"# emid: {(photometry['aper_pars']['tmid']+10*photometry['aper_errs']['dtt']):.6f} +- {photometry['aper_errs']['dt']:.6f}\n"
                # write to file with header
                output_file = os.path.join(obs, 'timeseries.csv')
                with open(output_file, 'w') as f:
                    f.write(file_header)
                    df.to_csv(f, index=False, float_format='%.6f')

                print(f"Saved {output_file}")

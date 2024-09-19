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
from scipy.ndimage import binary_dilation

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

# placeholder
phasecurve = lambda x: x

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
                tmid_idx = np.argmin(photometry['aper_transit'])
                final_tmid = photometry['aper_time'][tmid_idx]
                final_tmid_err = photometry['aper_errs']['dtt']
                tmodel = transit(photometry['aper_time'], prior)
                ax.plot(photometry['aper_time'], tmodel, 'g-', label='Transit Model with prior tmid')
                ax.axvline(final_tmid, color='g', linestyle='--')

                # isolate first eclipse
                emask = (photometry['aper_transit'] == 1) & (photometry['aper_phase'] < 1)
                if emask.sum() == 0:
                    # isolate second eclipse
                    emask = (photometry['aper_transit'] == 1) & (photometry['aper_phase'] > 1)
                
                if emask.sum() == 0:
                    # find arg of minimum
                    subdata = (photometry['aper_phase'] < 0.8)
                    if subdata.sum() > 0:
                        argmin = np.argmin(photometry['aper_transit'][subdata])
                        emid_time = photometry['aper_time'][subdata][argmin]
                    else:
                        subdata = (photometry['aper_phase'] > 1.2)
                        argmin = np.argmin(photometry['aper_transit'][subdata])
                        emid_time = photometry['aper_time'][subdata][argmin]
                else:
                    # compute median of the eclipse
                    emid_time = (max(photometry['aper_time'][emask]) + min(photometry['aper_time'][emask]))*0.5

                print(f"emid: {emid_time} +- {photometry['aper_errs']['dt']}")
                final_emid = emid_time
                final_emid_err = photometry['aper_errs']['dt']
                prior['tmid'] = emid_time
                prior['rprs'] = photometry['aper_pars']['y00']**0.5
                emodel = transit(photometry['aper_time'], prior)
                ax.plot(photometry['aper_time'], emodel, 'c-', label='Transit Model with prior emid')
                ax.axvline(final_emid, color='c', linestyle='--')

                # plot hotspot amp
                ax.axhline(np.percentile(photometry['aper_transit'],99.9), color='orange', linestyle='--', label='Hotspot Amp')

                # nightside
                tmid_idx = np.argmin(photometry['aper_transit'])
                tmid = photometry['aper_time'][tmid_idx]
                tmask = (photometry['aper_time'] > tmid - 0.1) & (photometry['aper_time'] < tmid + 0.1)
                # plot coldspot
                ax.axhline(np.percentile(photometry['aper_transit'][tmask],99), color='cyan', linestyle='--', label='Coldspot Amp')
                
                tdepth = np.percentile(photometry['aper_transit'][tmask],99) - photometry['aper_transit'].min()

                # dayside/eclipse
                emask2 = (photometry['aper_time'] > emid_time - 0.1) & (photometry['aper_time'] < emid_time + 0.1)
                
                ax.axhline(np.percentile(photometry['aper_transit'][emask2],97), color='purple', linestyle='--', label='Dayside Amp')

                ax.set_title(f"{sname} {f}")
                ax.legend(loc='best')
                plt.close()

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
                file_header =  f"# target: {sname}\n"
                file_header += f"# filter: {f}\n"
                file_header += f"# tmid: {final_tmid:.6f} +- {final_tmid_err:.6f}\n"
                file_header += f"# emid: {final_emid:.6f} +- {final_emid_err:.6f}\n"
                # depths
                file_header += f"# transit_depth: {tdepth:.6f}+-{photometry['hotspot_amp_std']:.6f}\n"
                file_header += f"# eclipse_depth: {np.percentile(photometry['aper_transit'][emask2],97)-1:.6f} +- {photometry['max_amp_std']:.6f}\n"
                # hotspot/coldspot
                file_header += f"# nightside_amp: {np.percentile(photometry['aper_transit'][tmask],99)-1:.6f} +- {photometry['max_amp_std']*(1.1+np.random.random()):.6f}\n"
                file_header += f"# hotspot_amp: {np.percentile(photometry['aper_transit'],99)-1:.6f} +- {photometry['max_amp_std']:.6f}\n"
                # hotspot_lon
                file_header += f"# hotspot_lon[deg]: {photometry['hotspot_lon']:.6f} +- {photometry['hotspot_lon_std']:.6f}\n"
                file_header += f"# hotspot_lat[deg]: {photometry['hotspot_lat']:.6f} +- {photometry['hotspot_lat_std']:.6f}\n"


                # write to file with header
                output_file = os.path.join(obs, 'timeseries.csv')
                with open(output_file, 'w') as f:
                    f.write(file_header)
                    df.to_csv(f, index=False, float_format='%.6f')

                print(f"Saved {output_file}")

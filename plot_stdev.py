import matplotlib
matplotlib.use('Agg')
import pickle
import glob
import shutil
import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats
from scipy.fft import fft, fftfreq

rsun = 6.955e8 # m
msun = 1.989e30 # kg
mjup = 1.898e27 # kg 
rjup = 7.1492e7 # m
mearth = 5.972e24 # kg
rearth = 6.3781e6 # m
au=1.496e11 # m 
G = 0.00029591220828559104 # day, AU, Msun

if __name__ == "__main__":
    phasecurve = lambda x: x
    targets = glob.glob("DATA/*")
    jet = plt.get_cmap('jet')

    print('computing stdev...')

    alldata = []
    datasets = []
    for t in targets:
        sname = t.lower().split('/')[-1]
        #with open('jwst_targets.txt') as f:
        #    jwst_targets = [line.rstrip().lower() for line in f]
        #if sname not in jwst_targets:
        #    continue

        filters = glob.glob(os.path.join(t,'IRAC*'))

        for f in filters:
            epochs = glob.glob(os.path.join(f,'*_phasecurve_0/'))
            # print(f,len(epochs))
            alltime = []
            allflux = []
            for e in epochs:
                
                    try:
                        prior = json.load(open(t+"/prior.json",'r'))
                        data = pickle.load(open(os.path.join(e,'data.pkl'),'rb'))
                    except:
                        continue

                    # photon noise calculation from in-transit data
                    tmask = data['aper_transit'] < 1

                    # convert MJy/sr to photons
                    try:
                        exptime = np.mean(data['exp_time'])
                    except:
                        exptime = np.diff(data['aper_time']).mean()*24*60*60 # seconds

                    gain = 3.7
                    flux_conv = 0.1257
                    photons = (data['aper_flux']/flux_conv)*gain*exptime
                    
                    # noise estimate
                    signal = np.sum(photons[tmask])
                    snr = np.sqrt(signal)
                    photon_noise_transit = 1/snr
                    photon_noise_timeseries = 1/np.sqrt(photons.mean())
                    
                    # transit depth error (rp/rs)^2
                    error_estimate = 2*data['aper_pars']['rprs']*data['aper_errs']['rprs']

                    noise_factor = error_estimate / photon_noise_transit
                    res_std = np.round(np.std(data['aper_residuals']/np.median(data['aper_flux'])),7)
                    nf_timeseries = res_std / photon_noise_timeseries

                    raw_residual = data['aper_flux']/np.median(data['aper_flux'])-data['aper_transit']
                    nf_timeseries_raw = np.std(raw_residual) / photon_noise_timeseries

                    dataset = {
                        'name':e.replace("DATA/","").split("/")[0],
                        'filter':f.split('/')[-1],
                        'cadence [sec]': np.round(np.diff(data['aper_time']).mean()*24*60*60,6),
                        'exp time [sec]':np.round(exptime,4),
                        'n observations':len(data['aper_time']),
                        'length of observations [day]': np.round(np.max(data['aper_time'])-min(data['aper_time']),4),
                        'photon_noise_factor_transit': np.round(noise_factor,4),
                        'photon_noise_factor_timeseries': np.round(nf_timeseries,4),
                        'photon_noise_factor_timeseries_raw': np.round(nf_timeseries_raw,4),
                        'average_photons_per_image':np.round(photons.mean(),1),
                        'residual_stdev': np.round(np.std(data['aper_residuals']/np.median(data['aper_flux'])),7),
                        'raw_stdev': np.round(np.std(raw_residual),7),
                    }

                    fig, ax = plt.subplots(3, figsize=(10,10))
                    nbins = int(np.sqrt(len(raw_residual)))
                    binspace = np.linspace(-0.02,0.02,201)
                    raw_label = f"Mean: {np.mean(raw_residual,):.4f} \n"\
                                f"Stdev: {np.std(raw_residual):.4f} \n"\
                                f"Skew: {stats.skew(raw_residual):.4f} \n"\
                                f"Kurtosis: {stats.kurtosis(raw_residual):.4f}\n"\
                                f"Photon Noise: {nf_timeseries_raw:.2f}"
                    ax[0].hist(raw_residual, bins=binspace,label=raw_label,color=jet(0.25),alpha=0.5)
                    rel_residuals = data['aper_residuals'] / np.median(data['aper_flux'])
                    detrend_label = f"Mean: {np.mean(rel_residuals):.4f} \n"\
                                f"Stdev: {np.std(rel_residuals):.4f} \n"\
                                f"Skew: {stats.skew(rel_residuals):.4f} \n"\
                                f"Kurtosis: {stats.kurtosis(rel_residuals):.4f}\n"\
                                f"Photon Noise: {nf_timeseries:.2f}"
                    ax[0].hist(rel_residuals, bins=binspace, label=detrend_label, color=jet(0.75),alpha=0.5)
                    ax[0].set_xlabel('Relative Flux Residuals')
                    ax[0].legend(loc='best')
                    ax[1].scatter(data['aper_time'], raw_residual, marker='.', label=f"Raw ({np.std(raw_residual,0)*100:.2f} %)",color=jet(0.25),alpha=0.25)
                    ax[1].scatter(data['aper_time'], rel_residuals, marker='.', label=f"Detrended ({np.std(rel_residuals,0)*100:.2f} %)",color=jet(0.75),alpha=0.25)
                    ax[1].legend(loc='best')
                    ax[1].set_xlabel('Time [BJD]')
                    ax[0].set_title(f'Residual Statistics: {e.replace("DATA/","")}')
                    ax[1].set_ylabel("Relative Flux")

                    # compute fourier transform of raw_residual
                    N = len(raw_residual)
                    fft_raw = fft(raw_residual)
                    fft_res = fft(rel_residuals)
                    xf = fftfreq(len(raw_residual), d=np.diff(data['aper_time']).mean()*24*60*60 )[:N//2]
                    #period = 1/xf
                    #mask = period < 0.1
                    fftraw = 2.0/N * np.abs(fft_raw[0:N//2])
                    #import pdb; pdb.set_trace()
                    # TODO sqaure + integrate under the curve and normalize such that it equals time series variance
                    # units of Hertz
                    #ax[2].semilogx(period*24, 2.0/N * np.abs(fft_raw[0:N//2]),alpha=0.5,label='Raw',color=jet(0.25))
                    #ax[2].semilogx(period*24, 2.0/N * np.abs(fft_res[0:N//2]),alpha=0.5,label='Detrended',color=jet(0.75))
                    ax[2].loglog(xf, 2.0/N * np.abs(fft_raw[0:N//2]),alpha=0.5,label='Raw',color=jet(0.25))
                    ax[2].loglog(xf, 2.0/N * np.abs(fft_res[0:N//2]),alpha=0.5,label='Detrended',color=jet(0.75))

                    ax[2].set_ylabel('Power')
                    #ax[2].set_xlabel('Period [hour]')
                    ax[2].set_xlabel('Frequency [Hz]')

                    #ax[2].set_xlim([0,2.5])
                    #ax[2].set_ylim([0,fftraw[mask].max()*1.05])
                    ax[2].legend()
                    ax[2].grid(True,ls='--')
                    plt.tight_layout()
                    plt.savefig(os.path.join(e,'residual_hist.png'))
                    plt.close()

                    # mask one hour of data
                    dt = np.diff(data['aper_time']).mean()
                    npts = int((10/24/60)/dt)
                    npts_1min = int((1/24/60)/dt)
                    ri = np.random.randint(data['aper_time'].shape[0])

                    std_flux = np.std(data['aper_flux'][ri:ri+npts])/np.median(data['aper_flux'][ri:ri+npts])
                    std_detrended = np.std(data['aper_detrended'][ri:ri+npts])/np.median(data['aper_detrended'][ri:ri+npts])

                    std_flux_1min = np.std(data['aper_flux'][ri:ri+npts_1min])/np.median(data['aper_flux'][ri:ri+npts_1min])

                    std_xcent = np.std(data['aper_xcent'][ri:ri+npts])#-np.median(data['aper_xcent'][ri:ri+npts])
                    std_xcent_1min = np.std(data['aper_xcent'][ri:ri+npts_1min])
                    std_ycent = np.std(data['aper_ycent'][ri:ri+npts])#-np.median(data['aper_ycent'][ri:ri+npts]))
                    std_ycent_1min = np.std(data['aper_ycent'][ri:ri+npts_1min])

                    fig, ax = plt.subplots(3,figsize=(9,10))
                    sub = np.min(data['aper_time'][ri:ri+npts]*24*60*60)
                    ax[0].scatter(data['aper_time'][ri:ri+npts]*24*60*60-sub, data['aper_flux'][ri:ri+npts]/np.median(data['aper_flux'][ri:ri+npts]), marker='.', label=f'Raw ({(100*std_flux):.2f} %)',color=jet(0.5),alpha=0.5)
                    ax[0].grid(True,ls='--')
                    #ax.set_title("Raw")
                    ax[0].scatter(data['aper_time'][ri:ri+npts]*24*60*60-sub, data['aper_detrended'][ri:ri+npts], marker='.', label=f'Detrended ({(100*std_detrended):.2f} %)',color=jet(0.75),alpha=0.5)
                    #ax.set_title("Detrended")
                    ax[0].set_xlabel('Time [sec]')
                    ax[0].set_ylabel('Relative Flux')
                    ax[0].set_title(f"10 minute sample: {e.replace('DATA/','')}")
                    # plot x and y centroid
                    ax[1].scatter(data['aper_time'][ri:ri+npts]*24*60*60-sub, data['aper_xcent'][ri:ri+npts], marker='.', label=f'X 5min ({std_xcent:.3f} px)',color=jet(0.25),alpha=0.5)
                    ax[1].scatter(data['aper_time'][ri:ri+npts]*24*60*60-sub, data['aper_xcent'][ri:ri+npts], marker='.', label=f'X 1min ({std_xcent_1min:.3f} px)',color=jet(0.25),alpha=0.5)
                    
                    ax[1].set_ylabel("X-Centroid [px]")
                    ax[1].set_xlabel('Time [sec]')
                    ax[2].set_xlabel('Time [sec]')
                    ax[2].scatter(data['aper_time'][ri:ri+npts]*24*60*60-sub, data['aper_ycent'][ri:ri+npts], marker='.', label=f'Y 5min ({std_ycent:.3f} px)',color=jet(0.75),alpha=0.5)
                    ax[2].scatter(data['aper_time'][ri:ri+npts]*24*60*60-sub, data['aper_ycent'][ri:ri+npts], marker='.', label=f'Y 1min ({std_ycent_1min:.3f} px)',color=jet(0.75),alpha=0.5)
                    
                    
                    ax[2].set_ylabel("Y-Centroid [px]")
                    ax[2].grid(True,ls='--')
                    ax[1].grid(True,ls='--')
                    ax[0].legend(loc='best')
                    ax[1].legend(loc='best')
                    ax[2].legend(loc='best')
                    plt.tight_layout()
                    plt.savefig(os.path.join(e,'10_minute_sample.png'))
                    plt.close()
                    
                    ri = np.random.randint(data['aper_time'].shape[0])
                    fig, ax = plt.subplots(3,figsize=(9,10))
                    sub = np.min(data['aper_time'][ri:ri+npts]*24*60*60)

#                     ax[0].scatter(1+data['aper_frames'][ri:ri+npts], data['aper_flux'][ri:ri+npts]/np.median(data['aper_flux'][ri:ri+npts]), marker='.', label=f'Raw',color=jet(0.5),alpha=0.5)
#                     ax[0].grid(True,ls='--')                        
#                     ax[0].set_xlabel('Frame Number')
#                     ax[0].set_ylabel('Relative Flux')
#                     ax[0].legend(loc='best')
#                     ax[0].set_title(f"frame_vs_centroid: {e.replace('DATA/','')}")
#                     # plot x and y centroid
#                     ax[1].scatter(1+data['aper_frames'][ri:ri+npts], data['aper_xcent'][ri:ri+npts], marker='.', label='X',color=jet(0.25),alpha=0.5)
#                     ax[1].set_ylabel("X-Centroid [px]")
#                     ax[1].set_xlabel('Frame Number')
#                     ax[2].set_xlabel('Frame Number')
#                     ax[2].scatter(1+data['aper_frames'][ri:ri+npts], data['aper_ycent'][ri:ri+npts], marker='.', label='Y',color=jet(0.75),alpha=0.5)
#                     for frame in np.unique(data['aper_frames']):
#                         fmask = data['aper_frames'][ri:ri+npts] == frame
#                         ax[1].scatter(1+data['aper_frames'][ri:ri+npts][fmask].mean(), data['aper_xcent'][ri:ri+npts][fmask].mean(), marker='o', label='X',color=jet(0.15),alpha=0.75)
#                         ax[2].scatter(1+data['aper_frames'][ri:ri+npts][fmask].mean(), data['aper_ycent'][ri:ri+npts][fmask].mean(), marker='o', label='Y',color=jet(0.85),alpha=0.75)
#                         ax[0].scatter(1+data['aper_frames'][ri:ri+npts][fmask].mean(), data['aper_flux'][ri:ri+npts][fmask].mean()/np.median(data['aper_flux'][ri:ri+npts]), marker='o', label='Raw',color=jet(0.35),alpha=0.95)
#                     ax[2].set_ylabel("Y-Centroid [px]")
#                     ax[2].grid(True,ls='--')
#                     ax[1].grid(True,ls='--')
#                     plt.tight_layout()
#                     plt.savefig(os.path.join(e,'frame_vs_centroid.png'))
#                     plt.close()
                    #import pdb; pdb.set_trace()
                    
                    fig, ax = plt.subplots(3, figsize=(10,10))
                    ax[0].scatter(data['aper_time'], data['aper_xcent'], marker='.',color=jet(0.25),alpha=0.5)
                    ax[0].set_xlabel('Time [BJD]')
                    ax[0].set_ylabel('X-Centroid [pixel]')
                    ax[1].scatter(data['aper_time'], data['aper_ycent'], marker='.',color=jet(0.75),alpha=0.5)
                    ax[1].set_xlabel('Time [BJD]')
                    ax[1].set_ylabel('Y-Centroid [pixel]')
                    # fft of x and y centroid
                    fft_x = fft(data['aper_xcent'])
                    fft_y = fft(data['aper_ycent'])
                    xf = fftfreq(len(data['aper_xcent']), d=np.diff(data['aper_time']).mean()*24*60*60 )[:N//2]
                    fft_x = 2.0/N * np.abs(fft_x[0:N//2])
                    fft_y = 2.0/N * np.abs(fft_y[0:N//2])
                    ax[2].loglog(xf, 2.0/N * np.abs(fft_x[0:N//2]),alpha=0.5,label='X-Centroid',color=jet(0.25))
                    ax[2].loglog(xf, 2.0/N * np.abs(fft_y[0:N//2]),alpha=0.5,label='Y-Centroid',color=jet(0.75))
                    ax[2].set_xlabel('Frequency [pixel/sec]')
                    ax[2].set_ylabel('Power')
                    ax[2].legend()
                    ax[2].grid(True,ls='--')
                    plt.tight_layout()
                    plt.savefig(os.path.join(e,'centroid_fft.png'))
                    plt.close()

                    json.dump(dataset,open(os.path.join(e,'residual_stats.json'),'w'),indent=4)
                    print(json.dumps(dataset,indent=4))
                    datasets.append(dataset)
    
    # make csv file from datasets
    df = pd.DataFrame(datasets)
    df.to_csv('residual_stats.csv',index=False)
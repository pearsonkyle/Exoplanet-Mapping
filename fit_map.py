import matplotlib
matplotlib.use('Agg')
import re
import os
import copy
import glob
import json
import pickle
import argparse
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from astropy import constants as const
from scipy import spatial
from astropy import units as u
from scipy import stats
from scipy.optimize import least_squares, brentq
from scipy.ndimage import gaussian_filter
import ultranest

from astropy.modeling.models import BlackBody
from scipy.interpolate import interp1d
from sklearn.feature_selection import mutual_info_regression
import starry
starry.config.lazy = False
starry.config.quiet = True

from tools import corner
from pylightcurve.models.exoplanet_lc import exotethys
from pylightcurve.models.exoplanet_lc import eclipse_mid_time

# TODO replace with astropy values
rsun = 6.955e8 # m
msun = 1.989e30 # kg
mjup = 1.898e27 # kg 
rjup = 7.1492e7 # m
mearth = 5.972e24 # kg
rearth = 6.3781e6 # m
au=1.496e11 # m 
pi = np.pi

def weightedflux(flux,gw,nearest):
    return np.sum(flux[nearest]*gw,axis=-1)

def gaussian_weights(X, w=None, neighbors=50, feature_scale=100):
    if isinstance(w, type(None)): w = np.ones(X.shape[1])
    Xm = (X - np.median(X,0))*w
    kdtree = spatial.cKDTree(Xm*feature_scale)
    nearest = np.zeros((X.shape[0],neighbors))
    gw = np.zeros((X.shape[0],neighbors),dtype=float)
    for point in range(X.shape[0]):
        ind = kdtree.query(kdtree.data[point],neighbors+1)[1][1:]
        dX = Xm[ind] - Xm[point]
        Xstd = np.std(dX,0)
        gX = np.exp(-dX**2/(2*Xstd**2))
        gwX = np.product(gX,1)
        gw[point,:] = gwX/gwX.sum()
        nearest[point,:] = ind
    gw[np.isnan(gw)] = 0.01
    gw[np.isinf(gw)] = 0.01
    return gw, nearest.astype(int)

# def phasecurve(time, prior):
#     # translate dictionary to starry api
# #     prior['system'].secondaries[0].map[1, :] = np.array([prior['y1m1'], prior['y10'], prior['y11']])
# #     prior['system'].secondaries[0].map.amp = prior['y00']
# #     prior['system'].secondaries[0].t0 = prior['tmid']
# #     prior['system'].secondaries[0].inc = prior['inc']
# #     return prior['system'].flux(time)
#     prior['planet'].map[1, :] = np.array([prior['y1m1'], prior['y10'], prior['y11']])
#     prior['planet'].map.amp = prior['y00']
#     prior['planet'].t0 = prior['tmid']
#     prior['planet'].inc = prior['inc']
#     # Instantiate the system
#     system = starry.kepler.System(prior['star'], prior['planet'])
#     flux = system.flux(time)
#     del system
#     return flux
# update fitting algorithm with different forward model calls
# try linear algebra solution


def mc_std(data,npts=100):
    # estimate stdev by cutting out a random sample of data, reduces outlier influence
    stds = []
    for i in range(100):
        stds.append(np.std(np.random.choice(data,npts)))
    return np.median(stds)

def spitzer_pixel_map(sv, title='', savedir=None):
    f,ax = plt.subplots(1,figsize=(8.5,7))
    im = ax.scatter(
        sv['aper_xcent'],
        sv['aper_ycent'],
        c=sv['aper_wf']/np.median(sv['aper_wf']),
        marker='.',
        vmin=0.97,
        vmax=1.03,
        alpha=0.5,
        cmap='jet',
    )
    # ax.set_xlim([
    #     np.median(sv['aper_xcent'])-3*mc_std(sv['aper_xcent']),
    #     np.median(sv['aper_xcent'])+3*mc_std(sv['aper_xcent'])
    # ])
    ax.set_xlim([
        np.percentile(sv['aper_xcent'],1)-0.1*mc_std(sv['aper_xcent']),
        np.percentile(sv['aper_xcent'],99)+0.1*mc_std(sv['aper_xcent'])
    ])
    ax.set_ylim([
        np.percentile(sv['aper_ycent'],1)-0.1*mc_std(sv['aper_ycent']),
        np.percentile(sv['aper_ycent'],99)+0.1*mc_std(sv['aper_ycent'])
    ])

    ax.set_title(title,fontsize=14)
    ax.set_xlabel('X-Centroid [px]',fontsize=14)
    ax.set_ylabel('Y-Centroid [px]',fontsize=14)
    cbar = f.colorbar(im)
    cbar.set_label('Relative Pixel Response',fontsize=14,rotation=270,labelpad=15)

    plt.tight_layout()
    if savedir:
        plt.savefig(savedir+title+".png")
        plt.close()
    else:
        return f

def spitzer_lightcurve(sv, savedir=None, suptitle=''):
    f,ax = plt.subplots(3,2,figsize=(12,12))
    f.suptitle(suptitle,y=0.99)
    res = sv['aper_residuals']/np.median(sv['aper_flux'])
    detrend = sv['aper_detrended']

    # #################### RAW FLUX ################
    ax[0,0].errorbar(
        sv['aper_time'], sv['aper_flux']/np.median(sv['aper_flux']),
        yerr=0,
        marker='.', ls='none', color='black',alpha=0.5)
    ax[0,0].set_xlim([min(sv['aper_time']), max(sv['aper_time'])])
    ax[0,0].set_xlabel('Time [JD]')
    ax[0,0].set_ylabel('Raw Relative Flux')
    ax[0,0].set_ylim([
        np.nanmean(detrend)-4*np.nanstd(detrend),
        np.nanmean(detrend)+4*np.nanstd(detrend)])

    # ################# DETRENDED FLUX ##################
    phase = ( sv['aper_time'] - sv['aper_pars']['tmid'] )/ sv['aper_pars']['per']
    ax[1,0].errorbar(
        phase, sv['aper_detrended'],
        yerr=0,
        marker='.', ls='none', color='black',alpha=0.15,
    )
    ax[1,0].plot(phase, sv['aper_transit'],'r-',zorder=4)
    bta,bfa = time_bin(sv['aper_time'], detrend)
    bpa = (bta - sv['aper_pars']['tmid'] )/ sv['aper_pars']['per']
    ax[1,0].plot(bpa, bfa, 'co', zorder=3, alpha=0.75)
    ax[1,0].set_xlim([min(phase), max(phase)])
    ax[1,0].set_xlabel('Phase')
    ax[1,0].set_ylabel('Relative Flux')
    ax[1,0].set_ylim([
        np.nanmean(detrend)-4*np.nanstd(detrend),
        np.nanmean(detrend)+4*np.nanstd(detrend)])

    # ################ RESIDUALS ###############
    bta,bfa = time_bin(sv['aper_time'], res)
    bstd = np.nanstd(bfa)*1e6
    std = np.nanstd(res)*1e6
    ax[2,0].plot(
        bta, bfa*1e6, 'co', zorder=3, alpha=0.75,
        label=r'$\sigma$ = {:.0f} ppm'.format(bstd)
    )
    ax[2,0].errorbar(
        sv['aper_time'], res*1e6,
        yerr=0,
        marker='.', ls='none', color='black',alpha=0.15,
        label=r'$\sigma$ = {:.0f} ppm'.format(std)
    )
    ax[2,0].legend(loc='best')
    ax[2,0].set_xlim([min(sv['aper_time']), max(sv['aper_time'])])
    ax[2,0].set_xlabel('Time [JD]')
    ax[2,0].set_ylabel('Residuals [ppm]')
    ax[2,0].set_ylim([
        np.nanmean(res*1e6)-3*np.nanstd(res*1e6),
        np.nanmean(res*1e6)+3*np.nanstd(res*1e6)])

    # ######## # # # # CENTROID X # # # #########
    ax[0,1].plot(
        sv['aper_time'], sv['aper_xcent'],
        marker='.', ls='none', color='black',alpha=0.5,
    )
    ax[0,1].set_xlim([min(sv['aper_time']), max(sv['aper_time'])])
    ax[0,1].set_xlabel('Time [JD]')
    ax[0,1].set_ylabel('X-Centroid [px]')
    ax[0,1].set_ylim([
        np.nanmean(sv['aper_xcent'])-3*np.nanstd(sv['aper_xcent']),
        np.nanmean(sv['aper_xcent'])+3*np.nanstd(sv['aper_xcent'])])

    ax[1,1].plot(
        sv['aper_time'], sv['aper_ycent'],
        marker='.', ls='none', color='black',alpha=0.5,
    )
    ax[1,1].set_xlim([min(sv['aper_time']), max(sv['aper_time'])])
    ax[1,1].set_xlabel('Time [JD]')
    ax[1,1].set_ylabel('Y-Centroid [px]')
    ax[1,1].set_ylim([
        np.nanmean(sv['aper_ycent'])-3*np.nanstd(sv['aper_ycent']),
        np.nanmean(sv['aper_ycent'])+3*np.nanstd(sv['aper_ycent'])])

    ax[2,1].plot(
        sv['aper_time'], sv['aper_npp'],
        marker='.', ls='none', color='black',alpha=0.5,
    )
    ax[2,1].set_xlim([min(sv['aper_time']), max(sv['aper_time'])])
    ax[2,1].set_xlabel('Time [JD]')
    ax[2,1].set_ylabel('Noise Pixel')
    ax[2,1].set_ylim([
        np.nanmean(sv['aper_npp'])-3*np.nanstd(sv['aper_npp']),
        np.nanmean(sv['aper_npp'])+3*np.nanstd(sv['aper_npp'])])

    plt.tight_layout(rect=[0, 0.03, 1, 0.98])
    if savedir:
        plt.savefig(savedir+suptitle+".png")
        plt.close()
    else:
        return f

def time_bin(time, flux, dt=2./(60*24)):
    bins = int(np.floor((max(time) - min(time))/dt))
    bflux = np.zeros(bins)
    btime = np.zeros(bins)
    for i in range(bins):
        mask = (time >= (min(time)+i*dt)) & (time < (min(time)+(i+1)*dt))
        if mask.sum() > 0:
            bflux[i] = np.nanmean(flux[mask])
            btime[i] = np.nanmean(time[mask])
    zmask = (bflux==0) | (btime==0) | np.isnan(bflux) | np.isnan(btime)
    return btime[~zmask], bflux[~zmask]

def sigma_clip(ogdata,dt, iterations=1):

    mask = np.ones(ogdata.shape,dtype=bool) & ~np.isnan(ogdata)
    for i in range(iterations):
        mdata = savgol_filter(ogdata[mask], dt, 2)
        res = ogdata[mask] - mdata
        try:
            std = np.nanmedian([np.nanstd(np.random.choice(res,100)) for i in range(250)])
        except:
            std = np.nanstd(res)
        mask[mask] = np.abs(res) < 3*std

    mdata = savgol_filter(ogdata[mask], dt, 2)
    res = ogdata[mask] - mdata

    data = deepcopy(ogdata)
    data[~mask] = np.nan
    return data, np.std(res)

def eclipse_ratio(priors,p='b',f='IRAC 3.6um'):
    # TODO wave as input?

    Te = priors['T*']*(1-0.1)**0.25 * np.sqrt(0.5/priors[p]['ars'])

    rprs = priors[p]['rp'] * const.R_jup / (priors['R*'] * const.R_sun)
    tdepth = rprs.value**2
    
    # bandpass integrated flux for planet
    wave36 = np.linspace(3.15,3.95,1000) * u.micron
    fplanet = BlackBody(Te*u.K)(wave36)
    fp36 = np.trapz(fplanet, wave36 )
    
    fstar = BlackBody(priors['T*']*u.K)(wave36)
    fs36 = np.trapz(fstar, wave36)

    wave45 = np.linspace(4,5,1000) * u.micron
    fplanet = BlackBody(Te*u.K)(wave45)
    fp45 = np.trapz(fplanet, wave45 )

    fstar = BlackBody(priors['T*']*u.K)(wave45)
    fs45 = np.trapz(fstar, wave45)
    
    print(" Stellar temp: {:.1f} K".format(priors['T*']))
    print(" Transit Depth: {:.4f} %".format(tdepth*100))
    print(" Equilibrium Temp: {:.0f} K".format(Te))
    if '3.6' in f:
        print(" Eclipse Depth @ IRAC 1 (3.6um): {:.0f} ppm".format(tdepth*fp36/fs36*1e6))
        print("         Fp/Fs @ IRAC 1 (3.6um): {:.4f}".format(fp36/fs36))
        return float(fp36/fs36)
    elif '4.5' in f:
        print(" Eclipse Depth @ IRAC 2 (4.5um): {:.0f} ppm".format(tdepth*fp45/fs45*1e6))
        print("         Fp/Fs @ IRAC 2 (4.5um): {:.4f}".format(fp45/fs45))
        return float(fp45/fs45)
    elif '5.8' in f:
        # bandpass integrated flux for planet
        wave = np.linspace(4.9,6.2,1000) * u.micron
        fplanet = BlackBody(Te*u.K)(wave)
        fp = np.trapz(fplanet, wave)
        fstar = BlackBody(priors['T*']*u.K)(wave)
        fs = np.trapz(fstar, wave)
        print(" Eclipse Depth @ IRAC 3 (5.8um): {:.0f} ppm".format(tdepth*fp/fs*1e6))
        print("         Fp/Fs @ IRAC 3 (5.8um): {:.4f}".format(fp/fs))
        return float(fp/fs)       
    elif '8.0' in f:
        # bandpass integrated flux for planet
        wave = np.linspace(6.5,9.1,1000) * u.micron
        fplanet = BlackBody(Te*u.K)(wave)
        fp = np.trapz(fplanet, wave)
        fstar = BlackBody(priors['T*']*u.K)(wave)
        fs = np.trapz(fstar, wave)
        print(" Eclipse Depth @ IRAC 4 (8.0um): {:.0f} ppm".format(tdepth*fp/fs*1e6))
        print("         Fp/Fs @ IRAC 4 (8.0um): {:.4f}".format(fp/fs))
        return float(fp/fs)

def plot_phasecurve(sv,return_temps=False):

    phase = (sv['aper_time']-sv['aper_pars']['tmid'])/sv['aper_pars']['per']
    bin_dt = 10./24./60.
    bt, bf = time_bin(sv['aper_time'], sv['aper_detrended'], bin_dt)
    bp = (bt-sv['aper_pars']['tmid'])/sv['aper_pars']['per']
    bt, br = time_bin(sv['aper_time'], sv['aper_residuals'], bin_dt)

    # Tb = brightnessTemp(sv['aper_pars'],sv['aper_filter'])
    # compute the brightness temperature over the orbit
    bcurve = brightness(bt,sv['aper_pars'])
    tbcurve = np.ones(bcurve.shape)
    for i,bc in enumerate(bcurve):
        sv['aper_pars']['fpfs'] = max((bc-1)/sv['aper_pars']['rprs']**2, 0.0001)
        try:
            tbcurve[i] = brightnessTemp(sv['aper_pars'],sv['aper_filter'])
        except:
            if i > 0:
                tbcurve[i] = tbcurve[i-1]
            continue

    if return_temps:
        return bt, bp, bf, tbcurve
        
    fig = plt.figure(figsize=(13,7))
    ax_lc = plt.subplot2grid((4,5), (0,0), colspan=5,rowspan=3)
    ax_res = plt.subplot2grid((4,5), (3,0), colspan=5, rowspan=1)
    axs = [ax_lc, ax_res]

    # residuals
    axs[1].plot(phase, sv['aper_residuals']/np.median(sv['aper_flux'])*1e6, 'k.', alpha=0.15, label=r'$\sigma$ = {:.0f} ppm'.format(np.std(sv['aper_residuals']/np.median(sv['aper_flux'])*1e6)))

    axs[1].plot(bp,1e6*br/np.median(sv['aper_flux']),'w.',zorder=2,label=r'$\sigma$ = {:.0f} ppm'.format(np.std(1e6*br/np.median(sv['aper_flux']))))

    axs[1].set_xlim([min(phase), max(phase)])
    axs[1].set_xlabel("Phase")
    axs[1].legend(loc='best')
    axs[1].set_ylabel("Residuals [ppm]")
    axs[1].grid(True,ls='--')

    axs[0].errorbar(phase, sv['aper_detrended'], yerr=np.std(sv['aper_residuals'])/np.median(sv['aper_flux']), ls='none', marker='.', color='black',alpha=0.1, zorder=1)

    # map color to equilibrium temperature
    im = axs[0].scatter(bp,bf,marker='o',c=tbcurve,vmin=500,vmax=3000,cmap='jet', zorder=2, s=20)
    cbar = plt.colorbar(im)
    cbar.ax.set_ylabel("B. Temp. [K]")

    axs[0].plot(phase, sv['aper_transit'], 'w--', zorder=3)
    axs[0].set_xlim([min(phase), max(phase)])
    axs[0].set_xlabel("Phase ")

    axs[0].set_ylabel("Relative Flux")
    axs[0].grid(True,ls='--')
    axs[0].set_ylim([0.965,1.025])

    plt.tight_layout()
    return fig

def brightnessTemp(priors,f='IRAC 3.6um'):
    # Solve for Tb using Fp/Fs, Ts and a filter bandpass
    if '3.6' in f or '36' in f:
        waveset = np.linspace(3.15, 3.9, 1000) * u.micron
    else:
        waveset = np.linspace(4,5,1000) * u.micron

    def f2min(T, *args):
        fpfs,tstar,waveset = args
        
        try:
            fplanet = BlackBody(T*u.K)(waveset)
            fstar = BlackBody(tstar*u.K)(waveset)
        except:
            fstar = BlackBody(waveset, tstar * u.K)
            fplanet = BlackBody(waveset, T * u.K)
        fp = np.trapz(fplanet, waveset)
        fs = np.trapz(fstar, waveset)
        return (fp/fs) - fpfs

    try:
        tb = brentq(f2min, 0.01,5500, args=(priors['fpfs'],priors['T*'],waveset))
    except:
        tb = -1
    return tb


class pc_fitter():

    def __init__(self, time, data, dataerr, prior, bounds, syspars, neighbors=100, mode='ns', verbose=False):
        self.time = time
        self.data = data
        self.dataerr = dataerr
        self.prior = copy.deepcopy(prior)
        self.bounds = bounds
        self.syspars = syspars
        self.neighbors = neighbors
        self.verbose = verbose
        
        if mode == 'ns':
            self.fit_nested()
        else:
            self.fit_lm()

    def fit_lm(self):

        freekeys = list(self.bounds.keys())
        boundarray = np.array([self.bounds[k] for k in freekeys])

        # trim data around predicted transit/eclipse time
        self.gw, self.nearest = gaussian_weights(self.syspars, neighbors=self.neighbors)

        def lc2min(pars):
            for i in range(len(pars)):
                self.prior[freekeys[i]] = pars[i]
            lightcurve = phasecurve(self.time, self.prior)
            detrended = self.data/lightcurve
            wf = weightedflux(detrended, self.gw, self.nearest)
            model = lightcurve*wf
            chi2 = ((self.data-model)/self.dataerr)**2
            return chi2

        res = least_squares(lc2min, x0=[self.prior[k] for k in freekeys], 
            bounds=[boundarray[:,0], boundarray[:,1]], jac='3-point', 
            loss='linear', verbose=True)

        self.parameters = copy.deepcopy(self.prior)
        self.errors = {}

        for i,k in enumerate(freekeys):
            self.parameters[k] = res.x[i]
            self.errors[k] = 0

        self.create_fit_variables()

    def fit_nested(self):
        freekeys = list(self.bounds.keys())
        boundarray = np.array([self.bounds[k] for k in freekeys])
        bounddiff = np.diff(boundarray,1).reshape(-1)

        # trim data around predicted transit/eclipse time
        self.gw, self.nearest = gaussian_weights(self.syspars, neighbors=self.neighbors)

        # local copy
        #system = starry.kepler.System(self.prior['star'], self.prior['planet'])

        def loglike(pars):
            for i in range(len(pars)):
                self.prior[freekeys[i]] = pars[i]
            lightcurve = self.prior['phasecurve'](self.time, self.prior)
            detrended = self.data/lightcurve
            wf = weightedflux(detrended, self.gw, self.nearest)
            model = lightcurve*wf

            residuals = self.data-model
            noise = np.std(residuals)/np.median(self.data)
            #fft = np.fft.fft(residuals)
            #N = len(residuals)
            #xf = np.fft.fftfreq(N, d=np.diff(self.time).mean()*24*60*60 )[:N//2]
            #power = 2.0/N * np.abs(fft[0:N//2])

            # ensure noise is greater than photon_noise
            penalty = len(residuals)*self.prior['photon_noise']**2 / noise**2
            chi2 = np.sum((residuals/self.dataerr)**2)

            #if (noise - self.prior['photon_noise']) < 0:
            #    import pdb; pdb.set_trace()

            return -0.5*(chi2)# + penalty)

        # create a gaussian priors, if any
        gaussian_priors = {}

        for k in ['y1m1']:
            if k in freekeys:
                gaussian_priors[k] = stats.norm(
                    self.bounds[k][0],
                    self.bounds[k][1])
                gaussian_priors[k+"_idx"] = freekeys.index(k)
            #def transform_1d(quantile):
            #    return gaussdistribution.ppf(quantile)

        def prior_transform(upars):
            vals = (boundarray[:,0] + bounddiff*upars)
            # set limits of phase amplitude to be less than eclipse depth
            # helps keep a positive flux on the night side
            if 'y00' in freekeys:
                edepth = vals[freekeys.index('y00')]
                for k in ['y10','y11']: # daynight, offset
                    if k in freekeys:
                        ki = freekeys.index(k)
                    else:
                        continue
                    if k == 'y10': # daynight (keep value positive)
                        vals[ki] = upars[ki]*edepth*0.75
                    #elif k == 'y11': # hot-spot offset
                    #    vals[ki] = upars[ki]*0.5*edepth # strictly positive or eastward
                        #upars[ki]*edepth*0.5-edepth*0.25 # [-0.25 -> 0.25] * eclipse depth
                        # could even constrain to be less than 50%?
                for k in gaussian_priors.keys():
                    if 'idx' in k:
                        continue
                    vals[gaussian_priors[k+"_idx"]] = gaussian_priors[k].ppf(upars[gaussian_priors[k+"_idx"]])
            return vals

        if self.verbose:
            self.results = ultranest.ReactiveNestedSampler(freekeys, loglike, prior_transform).run(max_ncalls=1e5)
        else:
            self.results = ultranest.ReactiveNestedSampler(freekeys, loglike, prior_transform).run(max_ncalls=1e5, show_status=self.verbose, viz_callback=self.verbose)

        self.errors = {}
        self.quantiles = {}
        self.parameters = copy.deepcopy(self.prior)

        for i, key in enumerate(freekeys):
            self.parameters[key] = self.results['maximum_likelihood']['point'][i]
            #if key == 'y11':
            #    #self.parameters[key] = self.results['posterior']['median'][i]
            #    diff = self.results['posterior']['errup'][i] - self.results['posterior']['errlo'][i]
            #    self.errors[key] = diff*0.5
            #else:
            self.errors[key] = self.results['posterior']['stdev'][i]
            self.quantiles[key] = [
                    self.results['posterior']['errlo'][i],
                    self.results['posterior']['errup'][i]]

        # self.results['maximum_likelihood']
        self.create_fit_variables()

    def create_fit_variables(self):
        self.transit = self.prior['phasecurve'](self.time, self.parameters)
        tmid = self.parameters['tmid'] - self.parameters['dt']
        self.phase = (self.time - tmid)/self.parameters['per']
        detrended = self.data / self.transit
        self.wf = weightedflux(detrended, self.gw, self.nearest)
        self.model = self.transit*self.wf
        self.detrended = self.data/self.wf
        self.detrendederr = self.dataerr 
        self.residuals = self.data - self.model
        self.chi2 = np.sum(self.residuals**2/self.dataerr**2)
        self.bic = len(self.bounds) * np.log(len(self.time)) - 2*np.log(self.chi2)

        # measured duration
        tdur = (self.transit < 1).sum() * np.median(np.diff(np.sort(self.time)))

    def plot_bestfit(self, bin_dt=10./(60*24), zoom=True, phase=True):
        f = plt.figure(figsize=(12,7))
        # f.subplots_adjust(top=0.94,bottom=0.08,left=0.07,right=0.96)
        ax_lc = plt.subplot2grid((4,5), (0,0), colspan=5,rowspan=3)
        ax_res = plt.subplot2grid((4,5), (3,0), colspan=5, rowspan=1)
        axs = [ax_lc, ax_res]

        bt, bf = time_bin(self.time, self.detrended, bin_dt)
        bp = (bt-self.parameters['tmid'])/self.parameters.get('per',self.prior['per'])

        if phase:
            axs[0].plot(bp,bf,'go',alpha=0.75,zorder=2)
            axs[0].plot(self.phase, self.transit, 'r-', zorder=3)
            axs[0].set_xlim([min(self.phase), max(self.phase)])
            axs[0].set_xlabel("Phase ")
        else:
            axs[0].plot(bt,bf,'go',alpha=0.75,zorder=2)
            axs[0].plot(self.time, self.transit, 'r-', zorder=3)
            axs[0].set_xlim([min(self.time), max(self.time)])
            axs[0].set_xlabel("Time [day]")
            
        axs[0].set_ylabel("Relative Flux")
        axs[0].grid(True,ls='--')

        if zoom:
            dmax = self.transit.max()-1
            axs[0].set_ylim([self.transit.min()-dmax, self.transit.max()+dmax])

        if phase:
            axs[0].errorbar(self.phase, self.detrended, yerr=np.std(self.residuals)/np.median(self.data), ls='none', marker='.', color='black', zorder=1, alpha=0.01)
        else:
            axs[0].errorbar(self.time, self.detrended, yerr=np.std(self.residuals)/np.median(self.data), ls='none', marker='.', color='black', zorder=1, alpha=0.01)

        bt, br = time_bin(self.time, self.residuals/np.median(self.data)*1e6, bin_dt)
        bp = (bt-self.parameters['tmid'])/self.parameters.get('per',self.prior['per'])

        if phase:
            axs[1].plot(self.phase, self.residuals/np.median(self.data)*1e6, 'k.', alpha=0.15, label=r'$\sigma$ = {:.0f} ppm'.format( np.std(self.residuals/np.median(self.data)*1e6)))
            axs[1].plot(bp,br,'g.',alpha=0.5,zorder=2,label=r'$\sigma$ = {:.0f} ppm'.format( np.std(br)))
            axs[1].set_xlim([min(self.phase), max(self.phase)])
            axs[1].set_xlabel("Phase")
        else:
            axs[1].plot(self.time, self.residuals/np.median(self.data)*1e6, 'k.', alpha=0.15, label=r'$\sigma$ = {:.0f} ppm'.format( np.std(self.residuals/np.median(self.data)*1e6)))
            axs[1].plot(bt,br,'g.',alpha=0.5,zorder=2,label=r'$\sigma$ = {:.0f} ppm'.format( np.std(br)))
            axs[1].set_xlim([min(self.time), max(self.time)])
            axs[1].set_xlabel("Time [day]")

        axs[1].legend(loc='best')
        axs[1].set_ylabel("Residuals [ppm]")
        axs[1].grid(True,ls='--')
        plt.tight_layout()
        return f,axs

    def plot_triangle(self):
        ranges = []
        mask1 = np.ones(len(self.results['weighted_samples']['logl']),dtype=bool)
        mask2 = np.ones(len(self.results['weighted_samples']['logl']),dtype=bool)
        mask3 = np.ones(len(self.results['weighted_samples']['logl']),dtype=bool)
        titles = []
        labels= []
        flabels = {
            'rprs':r'$\Delta$ (R$_{p}$/R$_{s}$)$^2$',
            'tmid':r'T$_{mid}$',
            'ars':r'a/R$_{s}$',
            'inc':r'Inc',
            'u1':r'u$_1$',
            'fpfs':r'F$_{p}$/F$_{s}$', 
            'omega':r'$\omega$',
            'ecc':r'$e$',
            'y00':'Uniform Brightness (y00)',
            'y10':'Day-Night (y10)',
            'y1m1':'Latitude (y1m1)',
            'y11':'Longitude (y11)',
            'dt':r'Mid-Eclipse $\Delta t$',
            'dtt':r'Mid-Transit $\Delta t$',
        }
        # constrain plots to +/- 4 sigma and estimate sigma levels
        for i, key in enumerate(self.quantiles):
            titles.append(f"{self.parameters[key]:.5f} +- {self.errors[key]:.5f}")

            if key == 'fpfs':
                ranges.append([
                    self.parameters[key] - 3*self.errors[key],
                    self.parameters[key] + 3*self.errors[key]
                ])
            elif key == 'y10': # day-night
                ranges.append([
                    max(0, self.parameters[key] - 3*self.errors[key]),
                    self.parameters[key] + 3*self.errors[key]
                ])
            else:
                ranges.append([
                    #self.parameters[key] - 4*self.errors[key],
                    #self.parameters[key] + 4*self.errors[key]
                    self.quantiles[key][0]-3*self.errors[key],
                    self.quantiles[key][1]+3*self.errors[key]
                ])

            mask3 = mask3 & (self.results['weighted_samples']['points'][:,i] > (self.parameters[key] - 3*self.errors[key]) ) & \
                (self.results['weighted_samples']['points'][:,i] < (self.parameters[key] + 3*self.errors[key]) )

            mask1 = mask1 & (self.results['weighted_samples']['points'][:,i] > (self.parameters[key] - self.errors[key]) ) & \
                (self.results['weighted_samples']['points'][:,i] < (self.parameters[key] + self.errors[key]) )

            mask2 = mask2 & (self.results['weighted_samples']['points'][:,i] > (self.parameters[key] - 2*self.errors[key]) ) & \
                (self.results['weighted_samples']['points'][:,i] < (self.parameters[key] + 2*self.errors[key]) )

            labels.append(flabels.get(key, key))

        chi2 = self.results['weighted_samples']['logl']*-2
        fig = corner(self.results['weighted_samples']['points'], 
            labels= labels,
            bins=int(np.sqrt(self.results['samples'].shape[0])), 
            range= ranges,
            #quantiles=(0.1, 0.84),
            plot_contours=False,
            levels=[chi2[mask1].max(), chi2[mask2].max(), chi2[mask3].max()],
            plot_density=False,
            titles=titles,
            data_kwargs={
                'c':chi2,
                'vmin':np.percentile(chi2[mask3],1),
                'vmax':np.percentile(chi2[mask3],99),
                'cmap':'viridis'
            },
            label_kwargs={
                'labelpad':42,
            },
            hist_kwargs={
                'color':'black',
            }
        )
        return fig


def get_ld(priors, band='irac1'):
    return exotethys(priors['LOGG*'], priors['T*'], priors['FEH*'], band, method='quad', stellar_model='atlas') 

def parse_args():
    parser = argparse.ArgumentParser()

    help_ = "Choose a target to process"
    parser.add_argument("-t", "--target", help=help_, type=str, default="HD189")

    help_ = "Choose a filter (3.6, 4.5, 5.8, 8.0)"
    parser.add_argument("-f", "--filter", help=help_, type=str, default="all")

    help_ = "Choose a planet (b, c, d, ...)"
    parser.add_argument("-p", "--planet", help=help_, type=str, default="all")

    help_ = "Fit only the eclipse"
    parser.add_argument("-e", "--eclipse", help=help_, default=False, action='store_true')

    help_ = "Fit only the transit"
    parser.add_argument("-tr", "--transit", help=help_, default=False, action='store_true')
    
    help_ = "Directory containing list of stars"
    parser.add_argument("-d", "--datadir", help=help_, default="DATA/", type=str)

    help_ = "Choose a frame number to process (0=all or 1-64)"
    parser.add_argument('-fr', '--frame', help=help_, default=0, type=int)

    help_ = "AOR observation filter"
    parser.add_argument("-a", "--aor", help=help_, default=0, type=int)

    return parser.parse_args()

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
            if args.target.lower() != os.path.basename(dirs[i][:-1].lower()):
                continue

        print(dirs[i])

        # TODO skip target if not present
        # load data
        try:
            priors = json.load(open(dirs[i]+"prior.json","r"))
            if args.aor == 0:
                photometry = pickle.load(open(dirs[i]+"photometry.pkl","rb"))
            else:
                photometry = pickle.load(open(dirs[i]+"photometry_"+str(args.aor)+".pkl","rb"))
        except:
            print(' no prior or no photometry files')
            continue
        # photometry{filter}.keys()
        # ['time', 'frame', 'aper_flux', 'aper_err', 'aper_xcent', 'aper_ycent', 'aper_npp', 'aper_bg', 
        # 'psf_flux', 'psf_err', 'psf_xcent', 'psf_ycent', 'psf_xsigma', 'psf_ysigma', 'psf_rot'])

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

        # alloc data
        sv = {}

        # for each filter
        for f in photometry.keys():

            #process specific filter
            if args.filter != 'all':
                if args.filter.lower() not in f:
                    continue

            # skip filter if no data
            if len(photometry[f]['time']) <= 1:
                print("skipping filter:",f)
                continue
            else:
                print("processing filter:",f)
            # convert data to numpy array
            for k in photometry[f]:
                photometry[f][k] = np.array(photometry[f][k])

            # loop through planets
            for p in priors['planets']:

                # https://arxiv.org/pdf/1001.2010.pdf eq 33
                w = priors[p].get('omega',0)
                tme = priors[p]['t0']+ priors[p]['period']*0.5 * (1 + priors[p]['ecc']*(4./np.pi)*np.cos(np.deg2rad(w)))
                tm = priors[p]['t0']

                period = priors[p]['period']
                offset = 0.25
                tphase = (photometry[f]['time'] - tm + offset*period)/period
                ephase = (photometry[f]['time'] - tme + offset*period)/period
                pdur = 2*np.arctan(priors['R*']*rsun/(priors[p]['sma']*au)) / (2*np.pi)
                # find transit events 
                events = np.unique(np.floor(tphase))
                for e in events:

                    tmask = (tphase > e + (offset-2*pdur)) & (tphase < e + (offset+2*pdur) )
                    mmask = (tphase > e + (offset-2*pdur+0.25)) & (tphase < e + (offset+2*pdur+0.25) )

                    emask = (ephase > e + (offset-1.25*pdur)) & (ephase < e + (offset+1.25*pdur) )
                    emask2 = (ephase > (e-1) + (offset-2*pdur)) & (ephase < (e-1) + (offset+2*pdur) ) # check for eclipse at prior orbit
                    tmask2 = (tphase > (e+1) + (offset-2*pdur)) & (tphase < (e+1) + (offset+2*pdur) ) # check for transit at next orbit
                    if args.eclipse:
                        if emask.sum() < 400:
                            print(f" not enough eclipse data ({emask.sum()}) at orbit #{e}")
                            continue

                        print("eclipse data detected")
                        print(" eclipse #{}: {}".format(e, emask.sum()))

                        mask = emask.copy()
                        tmid = priors[p]['t0'] + e*priors[p]['period']
                    elif args.transit:
                        if tmask.sum() < 400:
                            print(f" not enough transit data ({tmask.sum()}) at orbit #{e}")
                            continue

                        print("transit data detected")
                        print(" transit #{}: {}".format(e, tmask.sum()))

                        mask = tmask.copy()
                        tmid = priors[p]['t0'] + e*priors[p]['period']
                    else: # phasecurve
                        # skip if there's not enough data
                        if tmask.sum() < 400 or emask.sum() < 400 or mmask.sum() < 400:
                            if emask2.sum() > 400 and tmask.sum() > 400 and mmask.sum() > 400:
                                pass
                            else:
                                print(" not enough data ({},{},{}) at orbit #{}".format(tmask.sum(),mmask.sum(),emask.sum(),e))
                                continue
                        elif (tmask.sum() > 400 and emask.sum() > 400) or (tmask.sum() > 400 and emask2.sum() > 400):
                            pass
                        else:
                            print(" not enough data ({},{},{}) at orbit #{}".format(tmask.sum(),mmask.sum(),emask.sum(),e))
                            continue

                        print("phasecurve data detected")
                        print(" eclipse #{}: {}".format(e-1, emask2.sum()))
                        print(" transit #{}: {}".format(e, tmask.sum()))
                        print(" eclipse #{}: {}".format(e, emask.sum()))
                        print(" transit #{}: {}".format(e+1, tmask2.sum()))

                        mask = (tphase > e-1.75) & (tphase < e+1.25)

                        tmid = priors[p]['t0'] + e*priors[p]['period']
                        if emask2.sum() > 400: tmid = priors[p]['t0'] + (e-1)*priors[p]['period'] # adjust prior for transit model

                    # directory structure
                    # target/ filter/ midtransit
                    if args.eclipse:
                        savedir = os.path.join(dirs[i], f, str(int(tmid))+f"_eclipse_{args.frame}",)
                    elif args.transit:
                        savedir = os.path.join(dirs[i], f, str(int(tmid))+f"_transit_{args.frame}",)
                    else:
                        savedir = os.path.join(dirs[i], f, str(int(tmid))+f"_phasecurve_{args.frame}",)

                    if not os.path.exists( os.path.join(dirs[i], f) ):
                        os.mkdir( os.path.join(dirs[i],f) )

                    if not os.path.exists(savedir):
                        os.mkdir(savedir)

                    print(savedir)
                    # write aor to file
                    if args.aor == 0:
                        args.aor = min(photometry[f]['aor'][mask])

                    with open(os.path.join(savedir, f'aor_{args.aor}.txt'), 'w') as filename:
                        filename.write(str(args.aor))

                    time = photometry[f]['time'][mask]
                    smask = np.argsort(time)
                    dt = np.mean(np.diff(time[smask]))
                    ndt = max(11, int(5./(24*60*dt))*2+1)

                    # sigma clip flux before aperture selection
                    print('selecting aperture and sigma clipping...')
                    stds = []
                    
                    try:
                        for j in range(photometry[f]['aper_flux'].shape[1]):
                            # convert from photons to flux
                            photons = photometry[f]['aper_flux'][mask][smask][:,j]
                            exptime = photometry[f]['exptime'][mask][smask]
                            gain = photometry[f]['gain'][mask][smask]
                            flux_conv = photometry[f]['flux_conv'][mask][smask]
                            flux = photons/exptime/gain*flux_conv # MJy/sr
                            _, astd = sigma_clip(flux, ndt, iterations=2)
                            stds.append( astd)
                    except:
                        continue

                    # aperture selection
                    mi = np.argmin( stds )

                    # mask and sort photometry data
                    time = time[smask]
                    try:
                        photons, _ = sigma_clip( photometry[f]['aper_flux'][mask][smask][:,mi], ndt, iterations=2)
                        npp, _ = sigma_clip( photometry[f]['aper_npp'][mask][smask][:,mi], ndt, iterations=2)
                        wx, _ = sigma_clip( photometry[f]['aper_xcent'][mask][smask], ndt, iterations=2)
                        wy, _ = sigma_clip( photometry[f]['aper_ycent'][mask][smask], ndt, iterations=2)
                    except:
                        continue
                    frames = photometry[f]['frame'][mask][smask].astype(int)

                    #temp = photometry[f]['temp'][mask][smask]
                    # convert from photons to flux
                    exptime = photometry[f]['exptime'][mask][smask]
                    gain = photometry[f]['gain'][mask][smask]
                    flux_conv = photometry[f]['flux_conv'][mask][smask]
                    flux = photons/exptime/gain*flux_conv # MJy/sr
                    flux_err = np.sqrt(photons)/exptime/gain*flux_conv # MJy/sr

                    # remove all nan values + format for fitting
                    nanmask = np.isnan(time) | np.isnan(flux) | np.isnan(npp) | np.isnan(wx) | np.isnan(wy)
                    time = time[~nanmask]
                    flux = flux[~nanmask]
                    flux_err = flux_err[~nanmask]
                    npp = npp[~nanmask]
                    wx = wx[~nanmask]
                    wy = wy[~nanmask]
                    frames = frames[~nanmask]
                    exptime = exptime[~nanmask]

                    # remove first 30 min of data after any big gaps, ramp
                    tmask = np.ones(time.shape).astype(bool)

                    smask = np.argsort(time)
                    dts = np.diff(time[smask])
                    dmask = dts > (3/(24*60))
                    tmask[0:int(2*ndt)] = False # mask first 30 minutes of data
                    print("Gaps: ", dmask.sum())

                    idxs = np.argwhere(dmask).flatten()
                    if len(idxs) == 0:
                        ramp = np.exp(-np.arange(len(tmask))*dts[:-1].mean()/(0.69/24))
                        rramp = np.exp(-np.arange(len(tmask))[::-1]*dts[-1:].mean()/(0.69/24))
                    else:
                        ramp = np.exp(-np.arange(len(tmask))*dts[:idxs[0]].mean()/(0.69/24))
                        rramp = np.exp(-np.arange(len(tmask))[::-1]*dts[idxs[0]:].mean()/(0.69/24))

                        for idx in idxs:
                            ramp[idx-1:] = np.exp(-np.arange(len(ramp[idx-1:]))*dts[idx:].mean()/(0.69/24))

                        for idx in idxs[::-1]:  
                            rramp[:idx-1] = np.exp(-np.arange(len(rramp[:idx-1]))*dts[:idx].mean()/(0.69/24))[::-1]

                    #plt.plot(time, ramp, 'r-')
                    #plt.plot(time, rramp, 'b-')
                    #plt.savefig('test.png')
                    #plt.close()
                    #if not args.eclipse:
                    ramp += rramp

                    time = time[tmask]
                    flux = flux[tmask]
                    flux_err = flux_err[tmask]
                    npp = npp[tmask]
                    wx = wx[tmask]
                    wy = wy[tmask]
                    ramp = ramp[tmask] + np.random.random(tmask.sum())*0.01
                    frames = frames[tmask]
                    exptime = exptime[tmask]

                    # phase mask for WASP-121
                    if 'wasp-121' in sname:
                        phase = (time-tmid)/priors[p]['period']
                        if '3.6' in f:
                            # pmask = ((phase<0.65) | (phase>0.85))
                            pmask = (phase > 0.69) #& (phase < 1.2)
                        else:
                            pmask = np.ones(phase.shape[0],dtype=bool)

                        #pmask = np.zeros(phase.shape[0],dtype=bool)
                        #pmask[::10] = 1

                        time = time[pmask]
                        flux = flux[pmask]
                        flux_err = flux_err[pmask]
                        npp = npp[pmask]
                        wx = wx[pmask]
                        wy = wy[pmask]
                        ramp = ramp[pmask]
                        frames = frames[pmask]

                    # limb darkening from EXOFAST
                    if '3.6' in f:
                        lin,quad,_,_ = get_ld(priors,'irac1')
                    elif '4.5' in f:
                        lin,quad,_,_ = get_ld(priors,'irac2')
                    elif '5.8' in f:
                        lin,quad,_,_ = get_ld(priors,'irac3')
                    elif '8.0' in f:
                        lin,quad,_,_ = get_ld(priors,'irac4')

                    # define inputs
                    rprs = priors[p]['rp']*rjup/(priors['R*']*rsun)
                    inc = priors[p]['inc']
                    ars = priors[p].get('ars', priors[p]['sma']*au/(priors['R*']*rsun))

                    fpfs = eclipse_ratio(priors, p, f)

                    edepth = fpfs*1.25*rprs**2

                    # Instantiate the star
                    star = starry.Primary(starry.Map(ydeg=0, udeg=2, amp=1.0), m=priors['M*'], r=priors['R*'], prot=1.0)

                    star.map[1] = quad # quadratic
                    star.map[2] = lin # linear

                    # Instantiate the planet
                    planet = starry.kepler.Secondary(
                        starry.Map(ydeg=1, amp=edepth, reflected=False),  # the surface map
                        m=priors[p]['mass']*mjup/msun,  # mass in solar masses
                        r=priors[p]['rp']*rjup/rsun,  # radius in solar radii
                        inc=priors[p]['inc'],
                        porb=priors[p]['period'],  # orbital period in days
                        prot=priors[p]['period'],  # rotation period in days (synchronous)
                        ecc=priors[p]['ecc'],  # eccentricity
                        w=priors[p]['omega'],  # longitude of pericenter in degrees
                        t0=tmid,  # time of transit in days
                        theta0 = 180
                    )

                    tpars = {
                        'y00':edepth, # uniform brightness
                        'y10':0.33*edepth,     # day-night amplitude
                        'y11':0,         # hot-spot offset
                        'y1m1':0,        # latitude
                        'tmid':tmid,
                        'rprs':1, # scale factor for rprs
                        'dt':0,
                        'dtt':0, 
                        'per':priors[p]['period'],
                        'inc':priors[p]['inc'],
                    }

                    # only report params with bounds, all others will be fixed to initial value
                    mybounds = {
                        'rprs':[0,1.5], # factors of rprs
                        'dt':[-0.1,0.1],   # shift in mid-eclipse
                        'dtt':[-0.05, 0.05], # shift in mid-transit
                        'inc':[tpars['inc']-1.5, min(90,tpars['inc']+1.5)], # inclination
                        # spherical harmonics for map
                        'y00':[0.33*edepth,2*edepth],  # eclipse depth
                        'y10':[-edepth,edepth],  # day-night variation
                        'y1m1':[0,0.02*edepth],   # latitude - Gaussian Prior
                        'y11':[-0.1*edepth,0.2*edepth],    # hot-spot offset, positive is east
                    }

                    if args.eclipse:
                        del mybounds['rprs']
                        del mybounds['y10']
                        del mybounds['dtt']
                        #sdel mybounds['inc']
                    elif args.transit:
                        #mybounds['inc'] = [tpars['inc']-1.5, min(90,tpars['inc']+1.5)]
                        del mybounds['y10']
                        del mybounds['y00']
                        del mybounds['y11']
                        del mybounds['y1m1']
                        del mybounds['dt']
                    #else:
                    #    del mybounds['y1m1']
                    # if 'wasp-121' in sname:
                    #    phase = (time-tmid)/priors[p]['period']
                    #    if '3.6' in f:
                    #        del mybounds['y11']

                    # defaults to fitting a phasecurve for everything above

                    # Instantiate the system
                    system = starry.kepler.System(star, planet)

                    # HACK: run this to pre-compile the flux method
                    system.flux(0.0)

                    # upsample time array
                    super_time = np.linspace(time.min()-0.1, time.max()+0.1, 100000)

                    if 'inc' in mybounds:
                        # create inclination grids for spherical harmonics to interpolate over
                        incs = np.linspace(planet.inc-1,planet.inc+1,420)
                    else:
                        incs = np.linspace(planet.inc,planet.inc,1)

                    phasecurve_components = {
                        'tran_grid': np.zeros((len(incs), len(super_time))),
                        'y00_grid': np.zeros((len(incs), len(super_time))),
                        'y11_grid': np.zeros((len(incs), len(super_time))),
                        'y10_grid': np.zeros((len(incs), len(super_time))),
                        'y1m1_grid': np.zeros((len(incs), len(super_time)))
                    }

                    # precompute grid of transits at different inclinations
                    for jj,ii in enumerate(incs):
                        planet.inc = ii
                        
                        planet.map[1, :] = [0, 0, 0]
                        planet.map.amp = 0
                        phasecurve_components['tran_grid'][jj] = system.flux(super_time)

                        # basis set for eclipse map
                        planet.map[1, :] = [0, 0, 0]
                        planet.map.amp = 1
                        phasecurve_components['y00_grid'][jj]  = system.flux(super_time) - phasecurve_components['tran_grid'][jj] # uniform disk [flux scaled between 0 and 1]

                        # spherical harmonics for eclipse map
                        planet.map[1, :] = [1, 0, 0.0]
                        planet.map.amp = 1
                        phasecurve_components['y1m1_grid'][jj]  = system.flux(super_time) - phasecurve_components['y00_grid'][jj] - phasecurve_components['tran_grid'][jj] # latitude

                        planet.map[1, :] = [0, 1, 0.0]
                        planet.map.amp = 1
                        phasecurve_components['y10_grid'][jj]  = system.flux(super_time) - phasecurve_components['y00_grid'][jj] - phasecurve_components['tran_grid'][jj] # day-night

                        planet.map[1, :] = [0, 0, 1]
                        planet.map.amp = 1
                        phasecurve_components['y11_grid'][jj]  = system.flux(super_time) - phasecurve_components['y00_grid'][jj] - phasecurve_components['tran_grid'][jj] # longitude

                    planet.inc = priors[p]['inc']

                    if args.eclipse:
                        def phasecurve(times, pars):
                            mi = np.argmin(np.abs(incs-pars['inc']))
                            return 1+pars['y00']*np.interp(times+pars['dt'], super_time, phasecurve_components['y00_grid'][mi]) + \
                                pars['y1m1']*np.interp(times+pars['dt'], super_time, phasecurve_components['y1m1_grid'][mi]) + \
                                pars['y10']*np.interp(times+pars['dt'], super_time, phasecurve_components['y10_grid'][mi]) + \
                                pars['y11']*np.interp(times+pars['dt'], super_time, phasecurve_components['y11_grid'][mi])
                    elif args.transit:
                        def phasecurve(times, pars):
                            mi = np.argmin(np.abs(incs-pars['inc']))
                            return 1+pars['rprs']*np.interp(times+pars['dtt'], super_time, phasecurve_components['tran_grid'][mi]-1)
                    else:
                        def phasecurve(times, pars):
                            mi = np.argmin(np.abs(incs-pars['inc']))
                            # adding dt shifts model to the left
                            return 1+pars['rprs']*np.interp(times+pars['dtt'], super_time, phasecurve_components['tran_grid'][mi]-1)+\
                                pars['y00']*np.interp(times+pars['dt'], super_time, phasecurve_components['y00_grid'][mi]) + \
                                pars['y1m1']*np.interp(times+pars['dt'], super_time, phasecurve_components['y1m1_grid'][mi]) + \
                                pars['y10']*np.interp(times+pars['dt'], super_time, phasecurve_components['y10_grid'][mi]) + \
                                pars['y11']*np.interp(times+pars['dt'], super_time, phasecurve_components['y11_grid'][mi])

                    # forward model function
                    tpars['phasecurve'] = phasecurve

                    # mask by frame number in the datacube
                    if args.frame != 0:
                        fmask = frames >= args.frame
                        time = time[fmask] 
                        flux = flux[fmask] 
                        flux_err = flux_err[fmask]
                        exptime = exptime[fmask]
                        wx = wx[fmask]
                        wy = wy[fmask]
                        npp = npp[fmask]
                        ramp = ramp[fmask]

                    if args.eclipse:
                        #syspars = np.array([wx, wy, npp, ramp]).T
                        syspars = np.array([wx, wy, npp]).T
                    else:
                        syspars = np.array([wx, wy, npp, ramp]).T
                        #syspars = np.array([wx,wy,npp]).T

                    nneighbors = min(max(15,int(1/24/60/np.mean(np.diff(time)))),1000)

                    print(f"N Neighbors: {nneighbors}")
                    #import pdb; pdb.set_trace()
                    # make plot of priors
                    fig, ax = plt.subplots(3, figsize=(12,9))
                    ax[0].plot(time, flux, 'k.')
                    ax[0].set_xlabel("Time [BJD]")
                    ax[0].set_ylabel("Flux [MJy/sr]")
                    ax[1].plot(time, wx, 'r.',label='x-cent')
                    ax[1].plot(time, wy, 'b.',label='y-cent')
                    ax[1].set_ylabel('Pixel')
                    ax[1].set_xlabel("Time [BJD]")
                    tpars['dt'] = 0
                    model = phasecurve(time, tpars)
                    detrended = flux/model
                    gw, nearest = gaussian_weights(syspars, neighbors=nneighbors)
                    wf = weightedflux(detrended, gw, nearest)

                    tpars['dt'] = 0.1
                    model2 = phasecurve(time, tpars)
                    ax[2].plot(time, model2, 'b-') #should be shifted to the right

                    #ax[2].plot(time, flux/wf, 'k.',alpha=0.1)
                    btime, bflux = time_bin(time, flux/wf)
                    ax[2].plot(btime, bflux, 'go',alpha=0.9,label='Detrended + Binned')
                    ax[2].plot(time, model,'r-',label='Prior')
                    ax[2].legend(loc='best')
                    ax[1].legend(loc='best')
                    ax[2].set_ylim([model.min()-0.33*edepth, model.max()+0.33*edepth])
                    ax[2].set_xlabel("Time [BJD]")
                    ax[2].set_ylabel("Relative Flux")
                    plt.tight_layout()
                    plt.savefig(savedir+"/prior.png")
                    plt.close()

                    # photon noise estimate
                    # residuals should have more scatter than this
                    tpars['photon_noise'] = 1./np.sqrt(photons[~nanmask].mean())

                    # write data to file with columns: subt, aper, aper_err, syspars
                    #alldata = np.concatenate([time[:,None],flux[:,None],flux_err[:,None],syspars],axis=1)
                    #np.savetxt(f"wasp43b_transit.txt", alldata, header="time flux flux_err wxa wya npp ramp")
                    # write tpars to as json
                    #jdata = dict(**tpars)
                    #del jdata['phasecurve']
                    #with open(f"wasp43b_priors.json", 'w') as f: json.dump(jdata, f)
                    #import pdb; pdb.set_trace()


                    print("Running NS fit...")
                    myfit = pc_fitter(time, flux, flux_err, tpars, mybounds, syspars, mode='ns', verbose=True)

                    print('Values from NS fit:')
                    for k in myfit.bounds.keys():
                        print(" {}: {:.6f} +- {:.6f}".format(k, myfit.parameters[k], myfit.errors[k]))

                    rprs = planet.r/star.r
                    rprs2 = rprs**2                    
                    myfit.parameters['per'] = priors[p]['period']

                    #print("Plotting nsfit...")
                    #fig,axs = myfit.plot_bestfit(phase=True)
                    #plt.savefig(savedir+"/nsfit.png")
                    #plt.close()

                    myfit.plot_triangle()
                    plt.savefig(savedir+"/posterior.png")
                    plt.close()

                    # TODO improve data saving
                    #json.dump( {'parameters':sv['aper_pars'], 'errors':sv['aper_errs']}, open(savedir+"/parameters.json","w"),indent=4 )
                    #pickle.dump(sv,open(savedir+"/data.pkl","wb"))

                    #print("T$_{14}$ [BJD]  & Transit Duration   & %.5f $\pm$ %.5f \\\\" % ())

                    # make plot of bestfit
                    fig, ax = plt.subplots(5, figsize=(10,12))
                    ax[0].plot(time, flux, 'k.',alpha=0.5)
                    ndt = int(5./24/60/np.diff(time).mean())                    
                    sflux = gaussian_filter(flux,ndt)
                    ax[0].plot(time, sflux, 'w-',alpha=0.5)
                    ax[0].set_xlabel("Time [BJD]")
                    ax[0].set_ylabel("Flux [MJy/sr]")
                    ax[1].plot(time, wx, 'r.',label='x-cent',alpha=0.5)
                    ax[1].plot(time, wy, 'b.',label='y-cent',alpha=0.5)
                    ax[1].set_ylabel('Pixel')
                    ax[1].set_xlabel("Time [BJD]")
                    ax[2].plot(time, npp, 'k.',alpha=0.5)
                    #ax[2].plot(time, ramp, 'k-',alpha=0.5)
                    ax[2].set_ylabel('Noise Pixel')
                    ax[2].set_xlabel("Time [BJD]")

                    model = myfit.transit
                    detrended = flux/model
                    gw, nearest = gaussian_weights(syspars, neighbors=nneighbors)
                    wf = weightedflux(detrended, gw, nearest)
                    #ax[2].plot(time, flux/wf, 'k.',alpha=0.1)
                    if args.eclipse:
                        btime, bflux = time_bin(time, flux/wf)
                    else:
                        btime, bflux = time_bin(time, flux/wf, dt=10./(60*24))

                    ax[3].plot(btime, bflux, 'k.',alpha=0.9,label='Detrended + Binned')
                    tmid = myfit.parameters['dt'] + myfit.prior['tmid']
                    if args.eclipse:
                        ax[3].plot(time, model,'r-',label=f'E_mid = {tmid:.4f} +- {myfit.errors["dt"]:.4f} JD')
                    else:
                        ax[3].plot(time, model,'r-',label=f'Best fit')
                    # plot vertical line for emid
                    #ax[3].axvline(x=emid,ls='--',color='k')
                    ax[3].legend(loc='best')
                    ax[1].legend(loc='best')
                    ax[3].set_ylim([model.min()-0.33*edepth, model.max()+0.33*edepth])
                    ax[3].set_xlabel("Time [BJD]")
                    ax[3].set_ylabel("Relative Flux")
                    ax[3].grid(True,ls='--')
                    bmodel = phasecurve(btime, myfit.parameters)
                    ebar = 1e6*np.std(bflux-bmodel)
                    ax[4].errorbar(btime, 1e6*(bflux-bmodel),yerr=ebar,marker='.',ls='none',alpha=0.9,color='black', label=f"stdev = {ebar:.1f} ppm")
                    ax[4].axhline(y=0,ls='--',color='k')
                    ax[4].legend(loc='best')
                    ax[4].grid(True,ls='--')
                    ax[4].set_xlabel("Time [BJD]")
                    ax[4].set_ylabel("Residuals [ppm]")
                    plt.tight_layout()
                    plt.savefig(savedir+"/bestfit.png")
                    plt.close()


                    # make bestfit maps
                    planet.map[1, :] = [myfit.parameters['y1m1'],
                                        myfit.parameters['y10'],
                                        myfit.parameters['y11']]
                    planet.map.amp = myfit.parameters['y00']
                    if not args.transit:
                        # map of uncertainties on surface
                        maps = np.zeros((1500,1024,1024))
                        rmaps = np.zeros((1500,1024,1024))
                        fluxs = np.zeros(1500)
                        nfluxs = np.zeros(1500)
                        btemps = np.zeros(1500)
                        ntemps = np.zeros(1500)
                        mc_grids = {
                            'y00':[],
                            'y1m1':[],
                            'y10':[],
                            'y11':[],
                        }
                        tmask = (myfit.phase > 0.9) & (myfit.phase < 1.1)
                        for ri in range(1500):
                            y00 = np.random.normal(myfit.parameters['y00'], myfit.errors.get('y00',0))
                            y1m1 = np.random.normal(myfit.parameters['y1m1'], myfit.errors.get('y1m1',0))
                            y10 = np.random.normal(myfit.parameters['y10'], myfit.errors.get('y10',0))
                            y11 = np.random.normal(myfit.parameters['y11'], myfit.errors.get('y11',0))
                            scale = 1/y00
                            planet.map.reset()
                            planet.map[1, :] = [ 
                                scale*y1m1,
                                scale*y10,
                                scale*y11
                            ]
                            planet.map.amp = y00

                            dt = myfit.parameters['dt']
                            mi = np.argmin(np.abs(incs-myfit.prior['inc']))

                            mc_grids['y00'].append(planet.map.amp*np.interp(myfit.time+dt, super_time, phasecurve_components['y00_grid'][mi]))
                            mc_grids['y1m1'].append(y1m1*np.interp(myfit.time+dt, super_time, phasecurve_components['y1m1_grid'][mi]))
                            mc_grids['y10'].append(y10*np.interp(myfit.time+dt, super_time, phasecurve_components['y10_grid'][mi]))
                            mc_grids['y11'].append(y11*np.interp(myfit.time+dt, super_time, phasecurve_components['y11_grid'][mi]))
                            mc_flux = mc_grids['y00'][-1] + mc_grids['y1m1'][-1] + mc_grids['y10'][-1] + mc_grids['y11'][-1]
                            maps[ri] = planet.map.render(res=maps.shape[1])
                            rmaps[ri] = planet.map.render(projection='rect',res=maps.shape[1])

                            # dayside
                            flux = system.flux(time)
                            priors['fpfs'] = (flux.max() -1) / rprs**2
                            try:
                                btemps[ri] = brightnessTemp(priors, f=f)
                            except:
                                btemps[ri] = 0
    
                            fluxs[ri] = flux.max() # could use to compute brightness temperature

                            if not args.eclipse and not args.transit:
                                # nightside
                                nfluxs[ri] = np.nanmin(mc_flux[tmask])
                                priors['fpfs'] = (nfluxs[ri]) / rprs**2
                                if priors['fpfs'] < 0:
                                    ntemps[ri] = 0
                                else:
                                    ntemps[ri] = brightnessTemp(priors, f=f)

                        # maxfi = np.argmax(fluxs, axis=0)
                        # theta_max = theta[maxfi]
                        # flux_max = fluxs[:,maxfi]
                        # print(f"Hotspot at {-np.median(theta_max):.2f} +- {np.std(theta_max):.2f} degrees East")
                        # print(f"Hotspot flux at {np.median(flux_max)*1e6:.2f} +- {np.std(flux_max)*1e6:.2f} ppm")

                        # create distribution for long, lat of hotspot with likelihood at each point
                        lon_grid = np.linspace(-180, 180, maps[0].shape[0])
                        lat_grid = np.linspace(-90, 90, maps[0].shape[1])
                        lon_grid, lat_grid = np.meshgrid(lon_grid, lat_grid)
                        lon_grid = lon_grid.flatten()
                        lat_grid = lat_grid.flatten()

                        # find max of map within -90 to 90
                        max_lon = np.zeros(rmaps.shape[0])
                        max_lat = np.zeros(rmaps.shape[0])
                        max_amp = np.zeros(rmaps.shape[0])
                        min_amp = np.zeros(rmaps.shape[0])
                        mask = (lon_grid>-90) & (lon_grid<90) & (lat_grid>-90) & (lat_grid<90)
                        for ri in range(maps.shape[0]):
                            pflux = rmaps[ri].flatten()[mask]
                            mi = np.nanargmax(pflux)
                            max_lon[ri] = lon_grid[mask][mi]
                            max_lat[ri] = lat_grid[mask][mi]
                            max_amp[ri] = pflux[mi]
                            min_amp[ri] = np.nanmin(pflux)

                        # find max of mean
                        avg_map = np.median(rmaps,0)
                        max_i = np.nanargmax(avg_map.flatten()[mask])
                        max_lon_mean = lon_grid[mask][max_i]
                        max_lat_mean = lat_grid[mask][max_i]

                        omask = (max_lon>-80) & (max_lon<80) & (max_lat>-80) & (max_lat<80)

                        # histogram plot
                        fig, ax = plt.subplots(2,figsize=(5,8))
                        ax[0].hist(max_lon, bins=int(rmaps.shape[0]**0.5), range=(-45,45),label=f'Longitude: {max_lon_mean:.1f} +- {np.std(max_lon[omask]):.1f} deg',alpha=0.5,color='red')
                        ax[0].hist(max_lat, bins=int(rmaps.shape[0]**0.5), range=(-45,45),label=f'Latitude: {max_lat_mean:.1f} +- {np.std(max_lat[omask]):.1f} deg',alpha=0.5,color='blue')
                        ax[0].set_title("Planetary Hot Spot Location")
                        try:
                            ax[1].set_xlim([np.nanmedian(max_lon)-4*np.nanstd(max_lon[omask]),np.nanmedian(max_lon)+4*np.nanstd(max_lon[omask])])
                            ax[1].set_ylim([np.nanmedian(max_lat)-4*np.nanstd(max_lat[omask]),np.nanmedian(max_lat)+4*np.nanstd(max_lat[omask])])
                        except:
                            pass
                        ax[0].legend(loc='best')
                        ax[0].set_xlabel('Degrees')
                        ax[1].scatter(max_lon, max_lat, c='k',marker='.', alpha=0.25)
                        ax[1].set_xlabel('Longitude [Deg]')
                        ax[1].set_ylabel('Latitude [Deg]')
                        plt.tight_layout()
                        plt.savefig(savedir+"/max_lon_lat.png")
                        plt.close()

                        fig, ax = plt.subplots(1, figsize=(8, 4))
                        ax.set_xlim(-180, 180)
                        ax.set_ylim(-90, 90)
                        ax.grid(True,ls='--')
                        ax.set_xlabel('Longitude [Deg]')
                        ax.set_ylabel('Latitude [Deg]')
                        im = ax.imshow(np.median(rmaps,0), origin="lower", cmap="inferno", extent=(-180, 180, -90, 90), vmin=0,vmax=0.002)
                        cbar = plt.colorbar(im)
                        cbar.set_label("Intensity", rotation=270,labelpad=20, fontsize=14)
                        plt.tight_layout()
                        plt.savefig(savedir+"/bestmap_mean_rect.png")
                        plt.close()

                        fig, ax = plt.subplots(1, figsize=(8, 4))
                        ax.set_xlim(-180, 180)
                        ax.set_ylim(-90, 90)
                        ax.grid(True,ls='--')
                        ax.set_xlabel('Longitude [Deg]')
                        ax.set_ylabel('Latitude [Deg]')
                        im = ax.imshow(np.std(rmaps,0), origin="lower", cmap="inferno", extent=(-180, 180, -90, 90))
                        cbar = plt.colorbar(im)
                        cbar.set_label("Intensity Uncertainty", rotation=270,labelpad=20, fontsize=14)
                        plt.tight_layout()
                        plt.savefig(savedir+"/bestmap_std_rect.png")
                        plt.close()


                        fig, ax = plt.subplots(1, figsize=(7, 5))
                        ax.set_xlim(-1, 1)
                        ax.set_ylim(-1, 1)
                        ax.axis("off")
                        #edepth = 0.5*(myfit.transit.max()-1)
                        bestmap = np.median(maps,0)
                        im = ax.imshow(bestmap, origin="lower", cmap="inferno", extent=(-1, 1, -1, 1), vmin=0,vmax=0.002)
                        cbar = plt.colorbar(im)
                        cbar.set_label("Intensity", rotation=270,labelpad=20, fontsize=14)
                        plt.tight_layout()
                        plt.savefig(savedir+"/bestmap_mean.png")
                        plt.close()


                        # convert to temperature map
                        fbestmap = bestmap.flatten()
                        nanmask = ~np.isnan(fbestmap)
                        d2x = edepth / np.nansum(bestmap) # approximate integration factor
                        area = np.nansum(nanmask)
                        # Fp/Fs ~ d2x * sum(I)
                        uintensity = np.unique(fbestmap[nanmask])
                        relflux_range = np.linspace(uintensity.min(), uintensity.max(), 1000)*area*d2x/rprs2
                        temperature_range = np.zeros(relflux_range.shape)
                        for ii in range(len(relflux_range)):
                            priors['fpfs'] = relflux_range[ii] # estimate as uniform sphere of same intensity
                            temperature_range[ii] = brightnessTemp(priors)
                        ftemp = interp1d(relflux_range, temperature_range, kind='linear', bounds_error='False')
                        temp_map = ftemp(bestmap*area*d2x/rprs2)


                        # plot temperature map
                        fig, ax = plt.subplots(1, figsize=(7, 5))
                        ax.set_xlim(-1, 1)
                        ax.set_ylim(-1, 1)
                        ax.axis("off")
                        #edepth = 0.5*(myfit.transit.max()-1)
                        im = ax.imshow(temp_map, origin="lower", cmap="inferno", extent=(-1, 1, -1, 1), vmin=500,vmax=3000)
                        cbar = plt.colorbar(im)
                        cbar.set_label("Temperature [K]", rotation=270,labelpad=20, fontsize=14)
                        plt.tight_layout()
                        plt.savefig(savedir+"/tempmap_mean.png")
                        plt.close()

                        fig, ax = plt.subplots(1, figsize=(7, 5))
                        ax.set_xlim(-1, 1)
                        ax.set_ylim(-1, 1)
                        ax.axis("off")
                        beststd = np.std(maps,0)
                        im = ax.imshow(np.std(maps,0), origin="lower", cmap="inferno", extent=(-1, 1, -1, 1))
                        cbar = plt.colorbar(im)
                        cbar.set_label("Intensity Uncertainty", rotation=270,labelpad=20, fontsize=14)
                        plt.tight_layout()
                        plt.savefig(savedir+"/bestmap_std.png")
                        plt.close()

                        # convert to temperature map + uncertainty
                        fbestmap = (bestmap+beststd).flatten()
                        nanmask = ~np.isnan(fbestmap)
                        d2x = edepth / np.nansum(bestmap+beststd) # approximate integration factor
                        area = np.nansum(nanmask)
                        # Fp/Fs ~ d2x * sum(I)
                        uintensity = np.unique(fbestmap[nanmask])
                        relflux_range = np.linspace(uintensity.min(), uintensity.max(), 1000)*area*d2x/rprs2
                        temperature_range = np.zeros(relflux_range.shape)
                        for ii in range(len(relflux_range)):
                            priors['fpfs'] = relflux_range[ii] # estimate as uniform sphere of same intensity
                            temperature_range[ii] = brightnessTemp(priors)
                        ftemp = interp1d(relflux_range, temperature_range, kind='linear', bounds_error='False')
                        temp_map_up = ftemp((bestmap+beststd)*area*d2x/rprs2)

                        # convert to temperature map - uncertainty
                        fbestmap = (bestmap-beststd).flatten()
                        nanmask = ~np.isnan(fbestmap)
                        d2x = edepth / np.nansum(bestmap-beststd) # approximate integration factor
                        area = np.nansum(nanmask)
                        # Fp/Fs ~ d2x * sum(I)
                        uintensity = np.unique(fbestmap[nanmask])
                        relflux_range = np.linspace(uintensity.min(), uintensity.max(), 1000)*area*d2x/rprs2
                        temperature_range = np.zeros(relflux_range.shape)
                        for ii in range(len(relflux_range)):
                            priors['fpfs'] = relflux_range[ii] # estimate as uniform sphere of same intensity
                            try:
                                temperature_range[ii] = brightnessTemp(priors)
                            except:
                                pass
                        ftemp = interp1d(relflux_range, temperature_range, kind='linear', bounds_error='False')
                        temp_map_lo = ftemp((bestmap-beststd)*area*d2x/rprs2)

                        # plot temperature map
                        fig, ax = plt.subplots(1, figsize=(7, 5))
                        ax.set_xlim(-1, 1)
                        ax.set_ylim(-1, 1)
                        ax.axis("off")
                        #edepth = 0.5*(myfit.transit.max()-1)
                        im = ax.imshow(abs(temp_map_up-temp_map_lo), origin="lower", cmap="inferno", extent=(-1, 1, -1, 1),vmin=0,vmax=100)
                        cbar = plt.colorbar(im)
                        cbar.set_label("Temperature Uncertainty [K]", rotation=270,labelpad=20, fontsize=14)
                        plt.tight_layout()
                        plt.savefig(savedir+"/tempmap_std.png")
                        plt.close()


                    dt = myfit.parameters['dt']
                    mi = np.argmin(np.abs(incs-myfit.prior['inc']))

                    y00 = myfit.parameters['y00']*np.interp(myfit.time+dt, super_time, phasecurve_components['y00_grid'][mi])
                    y1m1 = myfit.parameters['y1m1']*np.interp(myfit.time+dt, super_time, phasecurve_components['y1m1_grid'][mi])
                    y10 = myfit.parameters['y10']*np.interp(myfit.time+dt, super_time, phasecurve_components['y10_grid'][mi])
                    y11 = myfit.parameters['y11']*np.interp(myfit.time+dt, super_time, phasecurve_components['y11_grid'][mi])

                    if not args.transit:
                        # create plot of each spherical harmonic
                        fig, ax = plt.subplots(2,2, figsize=(12, 9))
                        plt.subplots_adjust( left=0.075, right=0.95, top=0.95, bottom=0.05)
                        phase = myfit.phase%1
                        si = np.argsort(phase)
                        ax[0,0].plot(phase[si], y00[si]*1e6, 'k-', label='y00')
                        ax[0,0].set_xlabel("Phase")
                        ax[0,0].set_ylabel("Relative Flux [ppm]")
                        ax[0,0].set_title("y00 (Uniform Brightness)")
                        ax[0,0].grid(True,ls='--')
                        ax[0,1].plot(phase[si], y1m1[si]*1e6, 'k-', label='y1m1')
                        ax[0,1].set_title("y1m1 (Latitude)")
                        ax[0,1].set_xlabel("Phase")
                        ax[0,1].set_ylabel("Relative Flux [ppm]")
                        ax[0,1].grid(True,ls='--')
                        ax[1,0].plot(phase[si], y10[si]*1e6, 'k-', label='y10')
                        ax[1,0].set_title("y10 (Day-Night)")
                        ax[1,0].set_xlabel("Phase")
                        ax[1,0].set_ylabel("Relative Flux [ppm]")
                        ax[1,0].grid(True,ls='--')
                        ax[1,1].plot(phase[si], y11[si]*1e6, 'k-', label='y11')
                        ax[1,1].set_title("y11 (Longitude)")
                        ax[1,1].set_xlabel("Phase")
                        ax[1,1].set_ylabel("Relative Flux [ppm]")
                        ax[1,1].grid(True,ls='--')
                        y00_lower = np.percentile(mc_grids['y00'],16,axis=0)
                        y00_upper = np.percentile(mc_grids['y00'],100-16,axis=0)
                        y1m1_lower = np.percentile(mc_grids['y1m1'],16,axis=0)
                        y1m1_upper = np.percentile(mc_grids['y1m1'],100-16,axis=0)
                        y10_lower = np.percentile(mc_grids['y10'],16,axis=0)
                        y10_upper = np.percentile(mc_grids['y10'],100-16,axis=0)
                        y11_lower = np.percentile(mc_grids['y11'],16,axis=0)
                        y11_upper = np.percentile(mc_grids['y11'],100-16,axis=0)
                        ax[0,0].fill_between(phase[si], y00_lower[si]*1e6, y00_upper[si]*1e6, color='k', alpha=0.2)
                        ax[0,1].fill_between(phase[si], y1m1_lower[si]*1e6, y1m1_upper[si]*1e6, color='k', alpha=0.2)
                        ax[1,0].fill_between(phase[si], y10_lower[si]*1e6, y10_upper[si]*1e6, color='k', alpha=0.2)
                        ax[1,1].fill_between(phase[si], y11_lower[si]*1e6, y11_upper[si]*1e6, color='k', alpha=0.2)
                        fig.suptitle("Spherical Harmonic Components")
                        plt.savefig(savedir+"/bestfit_harmonics.png")
                        plt.close()

                        # write latex file to disk
                        # free parameters
                        sv['hotspot_lon'] = np.median(max_lon)
                        sv['hotspot_lat'] = np.median(max_lat)
                        sv['hotspot_lon_std'] = np.std(max_lon[omask])
                        sv['hotspot_lat_std'] = np.std(max_lat[omask])
                        sv['hotspot_amp'] = np.median(max_amp[omask])
                        sv['hotspot_amp_std'] = np.std(max_amp[omask])

                        sv['coldspot_amp'] = np.median(min_amp[omask])
                        sv['coldspot_amp_std'] = np.std(min_amp[omask])
                        sv['max_amp'] = np.median(fluxs)-1
                        sv['max_amp_std'] = np.std(fluxs)
                        sv['btemp'] = np.median(btemps)
                        sv['btemp_std'] = np.std(btemps)
                        sv['ntemp'] = np.median(ntemps)
                        sv['ntemp_std'] = np.std(ntemps)

                    # write data
                    sv['aper_time'] = myfit.time
                    sv['aper_flux'] = myfit.data # MJy/sr
                    sv['aper_err'] = myfit.dataerr
                    sv['aper_xcent'] = syspars[:,0]
                    sv['aper_ycent'] = syspars[:,1]
                    sv['aper_npp'] = syspars[:,2]
                    sv['aper_results'] = myfit.results
                    sv['aper_wf'] = myfit.wf
                    sv['aper_model'] = myfit.model
                    sv['aper_transit'] = myfit.transit
                    sv['aper_residuals'] = myfit.residuals
                    sv['aper_detrended'] = myfit.detrended
                    sv['aper_phase'] = myfit.phase
                    del myfit.parameters['phasecurve']
                    sv['aper_pars'] = myfit.parameters
                    sv['aper_errs'] = myfit.errors
                    sv['aper_frames'] = frames
                    sv['exp_time'] = exptime
                    sv['aor'] = args.aor
                    sv['lat_hist'] = max_lat
                    sv['lon_hist'] = max_lon
                    del tpars['phasecurve']
                    sv['aper_priors'] = tpars

                    try:
                        sv['aper_ramp'] = ramp
                    except:
                        pass

                    print("M$_s$ [M$_\odot$]   & Stellar Mass   & %.2f \\\\" % (star.m))
                    print("R$_s$ [R$_\odot$]   & Stellar Radius   & %.2f \\\\" % (star.r))
                    print("T$_s$ [K]   & Stellar Temperature  & %.2f \\\\" % (priors['T*']))
                    print("Fe/H   & Metallicity   & %.2f \\\\" % (priors['FEH*'])) 

                    print("\hline")
                    try:
                        rprs = planet.r/star.r
                        rprs2 = rprs**2 * myfit.parameters['rprs']
                        rprs2_err = rprs**2 * myfit.errors['rprs']
                        print("(R$_p$/R$_s$)$^2$   & Radius Ratio Squared   & %.5f $\pm$ %.5f \\\\" % (rprs2, rprs2_err))
                    except:
                        pass
                    print("Period [day]  & Orbital Period   & %.5f \\\\" % (myfit.prior['per']))

                    tmid = myfit.prior['tmid'] + myfit.parameters['dtt']
                    try:
                        print("T$_{mid}$ [BJD]  & Mid Transit Time   & %.5f $\pm$ %.5f \\\\" % (tmid, myfit.errors['dtt']))
                    except:
                        pass

                    #print("T$_{14}$ [BJD]  & Transit Duration   & %.5f $\pm$ %.5f \\\\" % ())
                    emid = eclipse_mid_time(myfit.prior['per'], 
                                            priors[p]['ars'], 
                                            priors[p]['ecc'], 
                                            myfit.parameters.get('inc', myfit.prior['inc']), 
                                            priors[p]['omega'], 
                                            myfit.prior['tmid']) + myfit.parameters['dt']
                    # does this need to be subtracted?!
                    sv['aper_emid'] = emid

                    try:
                        print("E$_{mid}$ [BJD]  & Mid Eclipse Time   & %.5f $\pm$ %.5f \\\\" % (emid, myfit.errors['dt']))
                    except:
                        pass
                    #print("E$_{14}$ [BJD]  & Transit Duration   & %.5f $\pm$ %.5f \\\\" % ())
                    try:
                        print("$i$ [deg]  & Inclination   & %.2f $\pm$ %.2f \\\\" % (myfit.parameters['inc'], myfit.errors['inc']))
                    except:
                        pass
                    print("$\omega$ [deg]  & Argument of Periapsis   & %.2f \\\\" % (planet.omega))
                    print("$e$ & Eccentricity   & %.2f \\\\" % (planet.ecc))

                    # spherical harmonic errors
                    print("Y$_00$ [ppm]   & Uniform Brightness   & %.2f $\pm$ %.2f \\\\" % (myfit.parameters['y00']*1e6, myfit.errors['y00']*1e6))
                    try:
                        print("Y$_10$ [ppm]   & Day-Night Amplitude   & %.2f $\pm$ %.2f \\\\" % (myfit.parameters['y10']*1e6, myfit.errors['y10']*1e6))
                    except:
                        pass
    
                    try:
                        print("Y$_11$ [ppm]   & Hot-spot Offset   & %.2f $\pm$ %.2f \\\\" % (myfit.parameters['y11']*1e6, myfit.errors['y11']*1e6))
                    except:
                        pass
                    try:
                        print("Y$_1m1$ [ppm]   & Hot-spot Latitude  & %.2f $\pm$ %.2f \\\\" % (myfit.parameters['y1m1']*1e6, myfit.errors['y1m1']*1e6))
                    except:
                        pass

                    # phase curve stuff
                    # hot-spot longitude
                    try:
                        print("&Hotspot Longitude [deg]   & %.2f $\pm$ %.2f \\\\" % (sv['hotspot_lon'], sv['hotspot_lon_std']))
                        print("&Hotspot Latitude [deg]   & %.2f $\pm$ %.2f \\\\" % (sv['hotspot_lat'], sv['hotspot_lat_std']))
                        print("&Eclipse Depth [ppm]   & %.2f $\pm$ %.2f \\\\" % (sv['max_amp']*1e6, sv['max_amp_std']*1e6))
                        print("&Day-side Temperature [K]   & %.0f $\pm$ %.0f \\\\" % (sv['btemp'], sv['btemp_std']))
                        print("&Night-side Temperature [K]   & %.0f $\pm$ %.0f \\\\" % (sv['ntemp'], sv['ntemp_std']))
                    except:
                        pass

                    # hot-spot amplitude
                    #print("Hotspot Intensity [ppm]   & %.2f $\pm$ %.2f \\\\" % (sv['hotspot_amp']*1e6, sv['hotspot_amp_std']*1e6))
                    # cold-spot amplitude
                    #print("Coldspot Intensity [ppm]   & %.2f $\pm$ %.2f \\\\" % (sv['coldspot_amp']*1e6, sv['coldspot_amp_std']*1e6))

                    resstd = sv['aper_residuals'].std() / sv['aper_flux'].mean()
                    print("&Residuals [ppm]   & %.2f  \\\\" % (resstd*1e6))
                    print("&Photon noise [ppm]   & %.2f \\\\" % (tpars['photon_noise']*1e6))
                    print("&Noise Factor  & %.2f \\\\" % (resstd/tpars['photon_noise']))
                    sv['noise_factor'] = resstd/tpars['photon_noise']
                    sv['photon_noise'] = tpars['photon_noise']

                    if args.frame == 0:
                        title = "/{} {}, {}, {}".format(
                            dirs[i].split('/')[1], p, f, "pixelmap"  )                 
                    else:
                        title = "/{} {}, {}, {}, {}".format(
                            dirs[i].split('/')[1], p, f, "pixelmap", args.frame  )                 

                    spitzer_pixel_map(sv, title, savedir=savedir)

                    if args.eclipse:
                        title = "/{} {}, {}, {}".format(dirs[i].split('/')[1], p, f, "eclipse" )
                    else:
                        title = "/{} {}, {}, {}".format(dirs[i].split('/')[1], p, f, "phasecurve" )
                    spitzer_lightcurve(sv, savedir=savedir, suptitle=title)

                    # see if residuals are still correlated to various parameters
                    fig,ax = plt.subplots(3, figsize=(4,7))
                    r = np.corrcoef(sv['aper_xcent'], sv['aper_residuals'])[0,1]
                    r2 = mutual_info_regression(myfit.syspars, sv['aper_residuals'])
                    ax[0].plot(sv['aper_xcent'], sv['aper_residuals']/sv['aper_flux'].mean(), 'k.',alpha=0.25, label=f"linear = {r:.3f}\n m.i. = {r2[0]:.3f}")
                    ax[0].set_xlabel("X Centroid")
                    ax[0].set_ylabel("Residuals")
                    ax[0].grid(True, ls='--')
                    ax[0].legend()

                    r = np.corrcoef(sv['aper_ycent'], sv['aper_residuals'])[0,1]
                    ax[1].plot(sv['aper_ycent'], sv['aper_residuals']/sv['aper_flux'].mean(),'k.',alpha=0.25, label=f"linear = {r:.3f}\n m.i. = {r2[1]:.3f}")
                    ax[1].set_xlabel("Y Centroid")
                    ax[1].set_ylabel("Residuals")
                    ax[1].grid(True, ls='--')
                    ax[1].legend()

                    r = np.corrcoef(sv['aper_npp'], sv['aper_residuals'])[0,1]
                    ax[2].plot(sv['aper_npp'], sv['aper_residuals']/sv['aper_flux'].mean(),'k.',alpha=0.25, label=f"linear = {r:.3f}\n m.i. = {r2[2]:.3f}")
                    ax[2].set_xlabel("Noise Pixel")
                    ax[2].set_ylabel("Residuals")
                    ax[2].grid(True, ls='--')
                    ax[2].legend()

                    fig.tight_layout()
                    fig.savefig(savedir+"/residuals.png", dpi=300)
                    plt.close()

                    sv['mutual_info_regression'] = r2

                    # pickle data
                    with open(savedir+"/data.pkl", 'wb') as pfile:
                        pickle.dump(sv, pfile)

                    # find residual correlations between npp, xcent and y cent with residuals

                    del myfit
                    try:
                        del maps, rmaps
                    except:
                        pass

                    #dict_keys(['hotspot_lon', 'hotspot_lat', 'hotspot_lon_std', 'hotspot_lat_std', 'hotspot_amp', 'hotspot_amp_std', 'coldspot_amp', 'coldspot_amp_std', 'max_amp', 'max_amp_std', 'btemp', 'btemp_std', 'ntemp', 'ntemp_std', 'aper_time', 'aper_flux', 'aper_err', 'aper_xcent', 'aper_ycent', 'aper_npp', 'aper_results', 'aper_wf', 'aper_model', 'aper_transit', 'aper_residuals', 'aper_detrended', 'aper_phase', 'aper_pars', 'aper_errs', 'aper_frames', 'exp_time', 'aor', 'lat_hist', 'lon_hist', 'aper_priors', 'aper_ramp', 'aper_emid', 'noise_factor', 'photon_noise', 'mutual_info_regression'])
                    # save a text file with time, detrended and error
                    # mask = sv['aper_phase'] < 0.6
                    # data = np.array([sv['aper_time'][mask], sv['aper_detrended'][mask], sv['aper_err'][mask]/sv['aper_flux'][mask]])
                    # np.savetxt("irac45_eclipse_1.txt", data.T, header="Time, Flux, Error")
                    # # plot
                    # from scipy.ndimage import uniform_filter1d
                    # plt.plot(sv['aper_time'][mask], sv['aper_detrended'][mask], 'k.')
                    # smooth = uniform_filter1d(sv['aper_detrended'][mask], size=1000)
                    # plt.plot(sv['aper_time'][mask], smooth, 'r-')
                    # plt.show()
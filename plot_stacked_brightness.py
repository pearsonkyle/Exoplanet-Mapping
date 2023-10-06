import pickle
import glob
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from ELCA_phasecurve import brightness, transit, eclipse

rsun = 6.955e8 # m
msun = 1.989e30 # kg
mjup = 1.898e27 # kg 
rjup = 7.1492e7 # m
mearth = 5.972e24 # kg
rearth = 6.3781e6 # m
au=1.496e11 # m 
G = 0.00029591220828559104 # day, AU, Msun

def phase_bin(time,flux,per,tmid=0,cadence=16,offset=0.25):
    '''
        Phase fold data and bin according to time cadence
        time - [days]
        flux - arbitrary unit
        per - period in [days]
        tmid - value in days
        cadence - spacing to bin data to [minutes] 
    '''
    phase = ((time-tmid)/per + offset)%1

    sortidx = np.argsort(phase)
    sortflux = flux[sortidx]
    sortphase = phase[sortidx]

    cad = cadence/60./24/per # phase cadence according to kepler cadence
    pbins = np.arange(0,1+cad,cad) # phase bins
    bindata = np.zeros(pbins.shape[0]-1)
    for i in range(pbins.shape[0]-1):
        pidx = (sortphase > pbins[i]) & (sortphase < pbins[i+1])

        if pidx.sum() == 0 or np.isnan(sortflux[pidx]).all():
            bindata[i] = np.nan
            continue

        bindata[i] = np.nanmedian(sortflux[pidx])

    phases = pbins[:-1]+np.diff(pbins)*0.5

    # remove nans
    #nonans = ~np.isnan(bindata)
    #return phases[nonans],bindata[nonans]
    return phases, bindata

def get_colors(inp, colormap, vmin=None, vmax=None):
    norm = plt.Normalize(vmin, vmax)
    return colormap(norm(inp))
    
if __name__ == "__main__":
    # function for loading certain pickles
    phasecurve = lambda x:x
    
    fig = plt.figure(figsize=(12,18))
    ax = [
        plt.subplot2grid((5, 2), (0, 0), rowspan=5), # 3.6
        plt.subplot2grid((5, 2), (0, 1), rowspan=5), # 4.5
    ]

    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes('right', size='5%', pad=0.05)

    offset = 0
    xoffset = 0
    targets = glob.glob("DATA/*")
    
    jet = plt.get_cmap('jet')
    cm = lambda x, minx,maxx: jet( (x-minx)/(maxx-minx) )
    
    print('computing temperatures...')
    temps = {}
    for t in targets:
        sname = t.lower().split('/')[-1]
        
        filters = glob.glob(os.path.join(t,'IRAC*'))

        for f in filters:
            epochs = glob.glob(os.path.join(f,'*_phasecurve_0'))
            # print(f,len(epochs))
            alltime = []
            allflux = []
            for e in epochs:
                prior = json.load(open(t+"/prior.json",'r'))
                temps[t] = prior['T*']*(1-0.3)**0.25 *(0.5*prior['b']['ars']**-1)**0.5
                try:
                    data = pickle.load(open(os.path.join(e,'data.pkl'),'rb'))
                except:
                    print(f'no data in {e}')
                    continue

    Tes = [temps[k] for k in temps]    
    sortidx = np.argsort(Tes)
    stargets = np.array(list(temps.keys()))[sortidx]
    for t in stargets:

        filters = glob.glob(os.path.join(t,'IRAC*'))
        for f in filters:
            epochs = glob.glob(os.path.join(f,'*_phasecurve_0'))
            print(f,len(epochs), np.round(temps[t]))
            if len(epochs) == 0:
                continue

            # acquire all data for global plot
            alltime = []
            allflux = []
            for e in epochs:
                prior = json.load(open(t+"/prior.json",'r'))
                try:
                    data = pickle.load(open(os.path.join(e,'data.pkl'),'rb'))
                except:
                    continue
                alltime.extend(data['aper_time'])
                allflux.extend(data['aper_detrended'])

            alltime = np.array(alltime)
            allflux = np.array(allflux)

            # plotting axis index
            if '3.6' in f:
                ai = 0
            elif '4.5' in f:
                ai = 1
            elif '5.8' in f:
                continue
                ai = 2
            elif '8.0' in f:
                continue
                ai = 3

            # plot raw data
            phase = (alltime - data['aper_pars']['tmid'] + 0.25* data['aper_pars']['per']) / data['aper_pars']['per'] % 1
            bphase, bflux = phase_bin(alltime, np.array(allflux), data['aper_pars']['per'], data['aper_pars']['tmid'], cadence=5)
            #ax[ai].plot(phase-0.25, allflux+offset, marker='.',color='gray', alpha=0.025, zorder =1)

            # global binning data
            bmask = (bflux > 1) & (bphase < 0.5)
            im = ax[ai].scatter((bphase-0.25)[bmask], (bflux+offset)[bmask], marker='.',cmap='jet',c=[temps[t]]*len(bflux[bmask]),vmin=500,vmax=3000, alpha=0.5, zorder=2)

            bmask = (bphase > 0.5)
            im = ax[ai].scatter((bphase-0.25)[bmask], (bflux+offset)[bmask], marker='.',cmap='jet',c=[temps[t]]*len(bflux[bmask]),vmin=500,vmax=3000, alpha=0.5, zorder=2)

            if ai == 1:
                ax[ai].text(-0.485,1.0+offset, t.split('/')[1]+" b", fontsize=9,color="black",zorder=4)
                #color=cm(temps[t],500,2500)
            
            # plot each epoch's phase curve
            for e in epochs:
                prior = json.load(open(t+"/prior.json",'r'))
                try:
                    data = pickle.load(open(os.path.join(e,'data.pkl'),'rb'))
                except:
                    continue
                phase = (data['aper_time'] - data['aper_pars']['tmid'] + 0.25* data['aper_pars']['per']) / data['aper_pars']['per'] % 1
                si = np.argsort(phase)

                # estimate transit duration and mast based on phase
                bright = data['aper_transit']
                pdur = 0.8*2*np.arctan(prior['R*']*rsun/(prior['b']['sma']*au)) / (2*np.pi)
                transit_mask = (phase-0.25 > -pdur) & (phase-0.25 < pdur)
                bright[transit_mask] = np.nan

                # plot phasecurve curve model
                ax[ai].plot(phase[si]-0.25, bright[si]+offset, ls='-', color='black', alpha=0.5, zorder=3)
                
        offset += 0.01

    cbar = fig.colorbar(im, cax=cax, orientation='vertical')
    cbar.ax.set_xlabel('Eq. Temp. [K]')
    ax[0].set_ylim([1-0.01,1+offset])
    ax[1].set_ylim([1-0.01,1+offset])
    ax[1].yaxis.set_ticklabels([])
    ax[0].set_xlabel("Phase")
    ax[1].set_xlabel("Phase")
    ax[0].set_ylabel("Relative Flux")
    ax[0].set_title("IRAC 3.6 um")
    ax[1].set_title("IRAC 4.5 um")
    ax[0].grid(True,ls='--')
    ax[1].grid(True,ls='--')
    plt.tight_layout()
    plt.subplots_adjust(wspace = 0.175)
    plt.savefig("FIGURES/brightness_curves.png")
    #plt.show()

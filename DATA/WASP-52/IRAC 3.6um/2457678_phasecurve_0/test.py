import pickle
from exotic.api.elca import transit
import json
import matplotlib.pyplot as plt

# Load the data
sv = pickle.load(open('data.pkl', 'rb'))

from pylightcurve import exotethys

priorf = json.load(open('../../prior.json', 'r'))

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

# dict_keys(['hotspot_lon', 'hotspot_lat', 'hotspot_lon_std', 'hotspot_lat_std', 'hotspot_amp', 'hotspot_amp_std', 'coldspot_amp', 'coldspot_amp_std', 'max_amp', 'max_amp_std', 'btemp', 'btemp_std', 'ntemp', 'ntemp_std', 'aper_time', 'aper_flux', 'aper_err', 'aper_xcent', 'aper_ycent', 'aper_npp', 'aper_results', 'aper_wf', 'aper_model', 'aper_transit', 'aper_residuals', 'aper_detrended', 'aper_phase', 'aper_pars', 'aper_errs', 'aper_frames', 'exp_time', 'aor', 'lat_hist', 'lon_hist', 'aper_priors', 'aper_ramp', 'aper_emid', 'noise_factor', 'photon_noise', 'mutual_info_regression'])

fig,ax = plt.subplots(1, figsize=(10,5))
ax.plot(sv['aper_time'], sv['aper_detrended'], 'k.', label='Detrended Flux')
ax.plot(sv['aper_time'], sv['aper_transit'], 'r-', label='Model')

# create a transit model
prior['tmid'] = sv['aper_pars']['tmid'] - sv['aper_pars']['dtt']
print(f"tmid: {prior['tmid']+prior['per']} +- {sv['aper_errs']['dtt']}")
tmodel = transit(sv['aper_time'], prior)
ax.plot(sv['aper_time'], tmodel, 'g-', label='Transit Model with prior tmid')
ax.axvline(prior['tmid']+prior['per'], color='g', linestyle='--')
prior['tmid'] = sv['aper_pars']['tmid'] +10*sv['aper_errs']['dtt'] + sv['aper_pars']['per']*0.5
print(f"emid: {prior['tmid']} +- {sv['aper_errs']['dt']}")
prior['rprs'] = sv['aper_pars']['y00']**0.5
emodel = transit(sv['aper_time'], prior)
ax.plot(sv['aper_time'], emodel, 'c-', label='Transit Model with prior emid')
ax.axvline(prior['tmid'], color='c', linestyle='--')
ax.legend(loc='best')
plt.show()

# write the timeseries to file
import pandas as pd 
df = pd.DataFrame({'time':sv['aper_time'], 'flux':sv['aper_detrended'], 'err': sv['aper_err']/sv['aper_flux'], 
                   'xcent':sv['aper_xcent'], 'ycent':sv['aper_ycent'], 'npp':sv['aper_npp'],
                   'phase':sv['aper_phase'], 'raw_flux':sv['aper_flux'], 'phasecurve':sv['aper_transit']})

# add tmid and emid to header with uncertainties
file_header = f"# tmid: {sv['aper_pars']['tmid']:.6f} +- {sv['aper_errs']['dtt']:.6f}\n"
file_header += f"# emid: {(sv['aper_pars']['tmid']+10*sv['aper_errs']['dtt']):.6f} +- {sv['aper_errs']['dt']:.6f}\n"
#df.to_csv('timeseries.csv', index=False)
# add header to the file
with open('timeseries.csv', 'w') as f:
    f.write(file_header)
    df.to_csv(f, index=False, float_format='%.6f')
import numpy as np
import astropy.units as u
from tqdm.auto import tqdm
import pandas as pd
import sys
import prospect.io.read_results as pread
from corner import quantile
from prospect.models import transforms
from prospect.sources import FastStepBasis
import matplotlib.pyplot as plt
import run_prosp_ml_sed
import fsps

sps = FastStepBasis()
#res, _ , _ = pread.results_from('galaxy_1_nonpara_fit.h5')
#sps = pread.get_sps(res)

sim_id = str(sys.argv[1])
snap = int(sys.argv[2])
gal_num = int(sys.argv[3])
full_phot = int(sys.argv[4])
z_idx = int(sys.argv[5])
#phot_file = str(sys.argv[4])
SNR = int(sys.argv[6])
try:
    extra_label = str(sys.argv[7])
    extra_label = f'_{extra_label}'
except:
    extra_label = ''
    print('no extra label')
    
if(full_phot == 1):
    prosp_out_folder = 'all_photometry'
elif(full_phot == 0):
    prosp_out_folder = 'limited_photometry'
else:
    raise Exception('not valid full_phot flag')

path_base = '/orange/narayanan/d.zimmerman/camels_results/'
phot_file = f'{path_base}/ml_data/all_filter_{sim_id}_snap{snap}.npz'
prosp_file = f'{path_base}/prosp_runs/{sim_id}/snap{snap}/{prosp_out_folder}/prosp_run_{sim_id}_snap{snap}_galaxy{gal_num}{extra_label}.h5'

outfile_prop = f'{path_base}/prosp_runs/{sim_id}/snap{snap}/{prosp_out_folder}/prosp_res_props_{sim_id}_snap{snap}_galaxy{gal_num}{extra_label}.npz'
outfile_spec = f'{path_base}/prosp_runs/{sim_id}/snap{snap}/{prosp_out_folder}/prosp_res_SEDs_{sim_id}_snap{snap}_galaxy{gal_num}{extra_label}.npz'

def get_sfh(res, mod):
    agebins = mod.params['agebins']
    thetas = mod.theta_labels()
    agebins_yrs = 10**agebins.T
    bin_edges = np.unique(agebins_yrs)
    dt = agebins_yrs[1, :] - agebins_yrs[0, :]
    epsilon = 1e-4 #fudge factor used to define the fraction time separation of adjacent points at the bin edges
    t = np.concatenate((bin_edges * (1.-epsilon), bin_edges * (1+epsilon)))
    t.sort()
    t = t[1:-1]
    zfrac_idx = [i for i, s in enumerate(thetas) if 'z_fraction' in s]
    zfrac_chain = res['chain'][:,zfrac_idx[0]:zfrac_idx[-1]+1]
    try:
        total_mass_chain = res['chain'][:,thetas.index('massmet_1')]
    except:
        total_mass_chain = res['chain'][:,thetas.index('logmass')]
    sfr_chain = []
    weights = res.get('weights',None)
    idx = np.argsort(weights)[-5000:]
    for i in idx:
        masses_chain = transforms.zfrac_to_masses(10**total_mass_chain[i], zfrac_chain[i], agebins)
        sfr = masses_chain / dt
        sfrout = np.zeros_like(t)
        sfrout[::2] = sfr
        sfrout[1::2] = sfr
        sfr_chain.append(sfrout)
    return (t[-1] - t[::-1])/1e9, sfr_chain


mass_quan = []
sfr_quan = []
metal_quan = []
dmass_quan = []
pd_phot, pd_sed=[],[]
phot_wave, pd_wave=[],[]
#res, obs, model = pread.results_from(f'/home/d.zimmerman/prospector_tutorial/g100_nonpara_backup.h5')
res, obs, model = pread.results_from(prosp_file)
print(model)
model = run_prosp_ml_sed.build_model(phot_file,z_idx)
print(model)
model_params = model.theta_labels()
spec, phot, mass_frac, dmass = [], [], [], []

#smass_true = np.load('const_norm_spec_gal0.npz',allow_pickle=True)['stellar_mass']
#Here I'm only calculating the quantiles of the 5000 most likely fit results

weights = res.get('weights',None)
idx = np.argsort(weights)[-5000:]


#if you want to process the entire posterior (which could really slow this calculation down, uncomment the next line
#idx = np.arange(len(res['chain']))
#print('-----')
#print(res['chain'][idx[0]])
#raise Exception()

for i in tqdm(idx):
    sspec, pphot, mmass_frac = model.predict(res['chain'][i], obs, sps)
    spec.append(sspec)
    phot.append(pphot)
    mass_frac.append(mmass_frac)
    dmass.append(sps.ssp.dust_mass)
raw_masses = res['chain'][idx,model_params.index('logmass')]
pred_masses = 10**raw_masses * mass_frac
mass_quan = quantile(np.log10(pred_masses), [.16, .5, .84], weights=res['weights'][idx])
dmass_quan = quantile(np.log10(dmass), [.16, .5, .84], weights=res['weights'][idx])
metal_quan = quantile(res['chain'][idx,model_params.index('logzsol')], [0.16, 0.5, 0.84], weights=res['weights'][idx])
raw_quan = quantile(raw_masses, [.16, .5, .84], weights=res['weights'][idx])
print()
print('stellar mass:',mass_quan)
print('dust mass:',dmass_quan)
print('metal val:',metal_quan)
print('mass formed:',quantile(raw_masses, [.16, .5, .84], weights=res['weights'][idx]))
print()
sfh_time, sfr_chain = get_sfh(res, model)
sfh_50, sfh_16, sfh_84 = [], [], []
#print(sfh_time)
#print(np.shape(sfh_time))
#print(sfr_chain)
#print(np.shape(sfr_chain))
#print(np.shape(res['weights']))
#print(idx)
for i in range(len(sfh_time)):
    sfh_quans = quantile([item[i] for item in sfr_chain], [.16, .5, .84], weights=res['weights'][idx])
    sfh_50.append(sfh_quans[1])
    sfh_16.append(sfh_quans[0])
    sfh_84.append(sfh_quans[2])

pd_phot.append(obs['maggies'])
phot_wave.append([x.wave_mean for x in obs['filters']])

pd_sed.append(obs['rest_sed'])
pd_wave.append(obs['wav'])

spec_50, spec_16, spec_84 = [], [], []
for i in range(len(spec[0])):
    quantiles = quantile([item[i] for item in spec], [.16, .5, .84], weights=res['weights'][idx])
    spec_50.append(quantiles[1])
    spec_16.append(quantiles[0])
    spec_84.append(quantiles[2])

"""
Now we save results to two pickle files, one with the SED info and one with the estimates for 
stellar mass, SFH, metallicity, and dust mass

The SED file contains the original Powderday SED, the photometry sampled from the Powderday SED (phot, phot_wave)
and the model Prospector SED

The masses and metallicities are saved as the median and the 16th-84th quantiles. same for the SFH and SED

"""



np.savez(outfile_prop,log_smass_quantiles = mass_quan, sfh_time = sfh_time, sfh_16 = sfh_16, sfh_50 = sfh_50, sfh_84 = sfh_84,
              logZsol = metal_quan, log_dmass_quantiles = dmass_quan,log_fmass_quantiles = raw_quan)
np.savez(outfile_spec,powderday_sed = pd_sed[0], powderday_wave = [item for item in pd_wave][0],
        spec_wave = sps.wavelengths*(1+np.load(phot_file)['z'][z_idx]),spec_50 = spec_50, spec_16 = spec_16, spec_84 = spec_84, phot = pd_phot, phot_wave = phot_wave)



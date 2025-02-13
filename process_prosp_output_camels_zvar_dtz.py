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
import run_prosp_ml_sed_zvar
import fsps

sps = FastStepBasis()
#res, _ , _ = pread.results_from('galaxy_1_nonpara_fit.h5')
#sps = pread.get_sps(res)
sim_id = str(sys.argv[1])
snap = int(sys.argv[2])
gal_num = int(sys.argv[3])
full_phot = int(sys.argv[4])
z_idx = int(sys.argv[5])
try:
    extra_label = str(sys.argv[6])
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

try:
    root_override = str(sys.argv[7])
    path_base = root_override
except:
    path_base = '/orange/narayanan/d.zimmerman/camels_results/'

phot_file = f'{path_base}/ml_data/all_filter_{sim_id}_snap{snap}.npz'
prosp_file = f'{path_base}/prosp_runs/{sim_id}/snap{snap}/{prosp_out_folder}/prosp_run_{sim_id}_snap{snap}_galaxy{gal_num}_z{z_idx}{extra_label}.h5'

outfile_prop = f'{path_base}/prosp_runs/{sim_id}/snap{snap}/{prosp_out_folder}/prosp_res_props_{sim_id}_snap{snap}_galaxy{gal_num}_z{z_idx}{extra_label}.npz'
outfile_spec = f'{path_base}/prosp_runs/{sim_id}/snap{snap}/{prosp_out_folder}/prosp_res_SEDs_{sim_id}_snap{snap}_galaxy{gal_num}_z{z_idx}{extra_label}.npz'

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

def get_sfh_z_var(res, mod):
    # pick out this stuff at the beginning
    weights = res.get('weights',None)
    idx = np.argsort(weights)[-5000:]
    thetas = mod.theta_labels()
    
    zfrac_idx = [i for i, s in enumerate(thetas) if 'z_fraction' in s]
    zfrac_chain = res['chain'][:,zfrac_idx[0]:zfrac_idx[-1]+1]
    try:
        total_mass_chain = res['chain'][:,thetas.index('massmet_1')]
    except:
        total_mass_chain = res['chain'][:,thetas.index('logmass')]

    time_chain = []
    sfr_chain = []
    for i in idx:
        agebins =  run_prosp_ml_sed_zvar.age_bin_calc([res['chain'][i,thetas.index('zred')]])
        agebins_yrs = 10**agebins.T
        bin_edges = np.unique(agebins_yrs)
        dt = agebins_yrs[1, :] - agebins_yrs[0, :]
        epsilon = 1e-4 #fudge factor used to define the fraction time separation of adjacent points at the bin edges
        t = np.concatenate((bin_edges * (1.-epsilon), bin_edges * (1+epsilon)))
        t.sort()
        t = t[1:-1]
        
        masses_chain = transforms.zfrac_to_masses(10**total_mass_chain[i], zfrac_chain[i], agebins)
        sfr = masses_chain / dt
        sfrout = np.zeros_like(t)
        sfrout[::2] = sfr
        sfrout[1::2] = sfr
        sfr_chain.append(sfrout)
        time_chain.append((t[-1] - t[::-1])/1e9)
    return time_chain, sfr_chain

mass_quan = []
sfr_quan = []
metal_quan = []
dmass_quan = []
pd_phot, pd_sed=[],[]
phot_wave, pd_wave=[],[]
#res, obs, model = pread.results_from(f'/home/d.zimmerman/prospector_tutorial/g100_nonpara_backup.h5')
res, obs, model = pread.results_from(prosp_file)
print(model)
model = run_prosp_ml_sed_zvar.build_model(phot_file,z_idx)
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
z_quan = quantile(res['chain'][idx,model_params.index('zred')], [0.16, 0.5, 0.84], weights=res['weights'][idx])
raw_quan = quantile(raw_masses, [.16, .5, .84], weights=res['weights'][idx])
z_red = res['chain'][idx,model_params.index('zred')]
print()
print('stellar mass:',mass_quan)
print('dust mass:',dmass_quan)
print('metal val:',metal_quan)
print('mass formed:',quantile(raw_masses, [.16, .5, .84], weights=res['weights'][idx]))
print('zred:',z_quan)
print('zmin:',np.min(z_red))
print('zmax:',np.max(z_red))
print()
sfh_time, sfr_chain = get_sfh_z_var(res, model)
print(np.shape(sfh_time))
print(np.shape(sfr_chain))
sfh_50, sfh_16, sfh_84 = [], [], []
#print(sfh_time)
#print(np.shape(sfh_time))
#print(sfr_chain)
#print(np.shape(sfr_chain))
#print(np.shape(res['weights']))
#print(idx)
print(sfh_time[:10])
for i in range(len(sfh_time[0])):
    sfh_quans = quantile([item[i] for item in sfr_chain], [.16, .5, .84], weights=res['weights'][idx])
    sfh_50.append(sfh_quans[1])
    sfh_16.append(sfh_quans[0])
    sfh_84.append(sfh_quans[2])

pd_phot.append(obs['maggies'])
phot_wave.append([x.wave_mean for x in obs['filters']])

pd_sed.append(obs['rest_sed'])
pd_wave.append(obs['wav'])

'''
for i in range(5000):
    plt.plot(sps.wavelengths*(1+z_red[i]),spec[i],alpha=0.1)
    #plt.axvline((sps.wavelengths*(1+z_red[i]))[0])#,color='black',linestyle='--',alpha=0.01)
plt.xscale('log')
plt.yscale('log')
plt.ylim(bottom=10**-20)
plt.savefig('/home/d.zimmerman/figures/test_zspec.png')
plt.close()
print(np.shape(spec))



for i in range(5000):
    plt.plot(sfh_time[i]-sfh_time[0][0],sfr_chain[i],alpha=0.1)
plt.savefig('/home/d.zimmerman/figures/test_zsfh.png')
plt.close()

diff_list = []
for i in range(5000):
    diff_list.append(sfh_time[i]-sfh_time[0])
plt.hist(np.ravel(diff_list),bins=np.arange(-0.105,0.115,0.01))
plt.savefig('/home/d.zimmerman/figures/test_zsfh_hist.png')
plt.close()

plt.hist(z_red)
plt.savefig('/home/d.zimmerman/figures/test_zsfh_z.png')
plt.close()
'''

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
              logZsol = metal_quan, log_dmass_quantiles = dmass_quan,log_fmass_quantiles = raw_quan,sfh_real = sfr_chain,z_real = z_red,z_quan = z_quan)
np.savez(outfile_spec,powderday_sed = pd_sed[0], powderday_wave = [item for item in pd_wave][0],
        spec_wave = sps.wavelengths,spec_50 = spec_50, spec_16 = spec_16, spec_84 = spec_84, phot = pd_phot, phot_wave = phot_wave,spec_real = )


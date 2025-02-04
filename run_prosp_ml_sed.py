import numpy as np
from sedpy.observate import load_filters
import h5py
import prospect.io.read_results as pread
from prospect.models import priors, transforms,sedmodel
from prospect.sources import FastStepBasis
from scipy.stats import truncnorm
from prospect.io import write_results as writer
from prospect.fitting import fit_model
import sys, os
#from astropy.cosmology import Planck13
from astropy.cosmology import FlatLambdaCDM
from hyperion.model import ModelOutput
from astropy import units as u
from astropy import constants


#------------------------
# Convienence Functions
#------------------------

def get_best(res, **kwargs):
    imax = np.argmax(res['lnprobability'])
    theta_best = res['chain'][imax, :].copy()
    return theta_best

def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx

def zfrac_to_masses_log(logmass=None, z_fraction=None, agebins=None, **extras):
    sfr_fraction = np.zeros(len(z_fraction) + 1)
    sfr_fraction[0] = 1.0 - z_fraction[0]
    for i in range(1, len(z_fraction)):
        sfr_fraction[i] = np.prod(z_fraction[:i]) * (1.0 - z_fraction[i])
    sfr_fraction[-1] = 1 - np.sum(sfr_fraction[:-1])
    # convert to mass fractions
    time_per_bin = np.diff(10**agebins, axis=-1)[:, 0]
    mass_fraction = sfr_fraction * np.array(time_per_bin)
    mass_fraction /= mass_fraction.sum()
    
    if (mass_fraction < 0).any():
        idx = mass_fraction < 0
        if np.isclose(mass_fraction[idx],0,rtol=1e-8):
            mass_fraction[idx] = 0.0
        else:
            raise ValueError('The input z_fractions are returning negative masses!')
        
    masses = 10**logmass * mass_fraction
    return masses

#----------------------
# SSP function
#-----------------------

def build_sps(**kwargs):
    """
    This is our stellar population model which generates the spectra for stars of a given age and mass.
    Because we are using a non parametric SFH model, we do have to use a different SPS model than before
    """
    from prospect.sources import FastStepBasis
    sps = FastStepBasis(zcontinuous=1)
    return sps



def build_model(spec_dir,z_idx,**kwargs):

    z = np.load(spec_dir)['z'][z_idx]
    print('z:',z)
    cosmo = FlatLambdaCDM(Om0=0.3,Tcmb0 = 2.725,H0 = 67.11)
    if(z<10**-4):
        dl = (10*u.Mpc)
    else:
        dl = cosmo.luminosity_distance(z).to(u.Mpc)
    
    print('building model')
    model_params = []
    #basics
    model_params.append({'name': "lumdist", "N": 1, "isfree": False,"init": dl.value,"units": "Mpc"})
    model_params.append({'name':'zred','N':1,'isfree':False,'init':z})

    model_params.append({'name': 'imf_type', 'N': 1,'isfree': False,'init': 1})
    model_params.append({'name': 'dust_type', 'N': 1,'isfree': False,'init': 2,'prior': None})
    model_params.append({'name': 'dust2', 'N': 1,'isfree': True, 'init': 0.1,'prior': priors.ClippedNormal(mini=0.0, maxi=2.0, mean=0.0, sigma=0.3)})
    #model_params.append({'name': 'dust2', 'N': 1,'isfree': False, 'init': 0.0,'prior': None})
    model_params.append({'name': 'add_dust_emission', 'N': 1,'isfree': False,'init': 1,'prior': None})
    model_params.append({'name': 'duste_gamma', 'N': 1,'isfree': True,'init': 0.01,'prior': priors.TopHat(mini=0.0, maxi=1.0)})
    model_params.append({'name': 'duste_umin', 'N': 1,'isfree': True,'init': 1.0,'prior': priors.TopHat(mini=0.1, maxi=20.0)})
    model_params.append({'name': 'duste_qpah', 'N': 1,'isfree': True,'init': 3.0,'prior': priors.TopHat(mini=0.0, maxi=6.0)})
    
    model_params.append({'name': 'add_agb_dust_model', 'N': 1,'isfree': False,'init': 0})
    
    #M-Z
    model_params.append({'name': 'logmass', 'N': 1,'isfree': True,'init': 8.0,'prior': priors.Uniform(mini=7., maxi=12.)})
    
    # CHANGE
    model_params.append({'name': 'logzsol', 'N': 1,'isfree': True,'init': -0.5,'prior': priors.Uniform(mini=-1.5, maxi=0.5)})
    #model_params.append({'name': 'logzsol', 'N': 1,'isfree': False,'init': 0,'prior': None})

    model_params.append({'name': "sfh", "N": 1, "isfree": False, "init": 3})
    model_params.append({'name': "mass", 'N': 3, 'isfree': False, 'init': 1., 'depends_on':zfrac_to_masses_log})
    model_params.append({'name': "agebins", 'N': 1, 'isfree': False,'init': []})
    model_params.append({'name': "z_fraction", "N": 2, 'isfree': True, 'init': [0, 0],'prior': priors.Beta(alpha=1.0, beta=1.0, mini=0.0, maxi=1.0)})       
    
    
    #here we set the number and location of the timebins, and edit the other SFH parameters to match in size
    n = [p['name'] for p in model_params]
    tuniv = np.round(cosmo.age(z).to('Gyr').value,decimals=1) #Gyr, age at z=0 #CHANGE
    #space = 0.1
    
    agelims = np.linspace(8,np.log10(tuniv)+9,9)
    agelims = [0]+agelims.tolist()
    
    print(agelims)
    nbins = len(agelims)-1
    agebins = np.array([agelims[:-1], agelims[1:]])
    #print(agebins)
    #raise Exception()
    zinit = np.array([(i-1)/float(i) for i in range(nbins, 1, -1)])
    #zinit = np.array([0.1 for i in range(nbins, 1, -1)])
    
    #plt.plot(
    # Set up the prior in `z` variables that corresponds to a dirichlet in sfr
    # fraction.
    alpha = np.arange(nbins-1, 0, -1)
    zprior = priors.Beta(alpha=alpha, beta=np.ones_like(alpha), mini=0.0, maxi=1.0)
    
    model_params[n.index('mass')]['N'] = nbins
    model_params[n.index('agebins')]['N'] = nbins
    model_params[n.index('agebins')]['init'] = agebins.T
    model_params[n.index('z_fraction')]['N'] = nbins-1
    model_params[n.index('z_fraction')]['init'] = zinit
    model_params[n.index('z_fraction')]['prior'] = zprior
    
    model = sedmodel.SedModel(model_params)
    
    return model


#------------------
# Build Observations
#-------------------

def build_obs(spec_dir,gnum,sn,z_idx,filt_file,**kwargs):
    obs = {}
    dat = np.load(spec_dir)
    all_filters = dat['filter_names']
    
    filt_list = []
    filt_mask = []
    
    if(filt_file == 'n/a'):
        filt_mask = np.arange(0,len(all_filters),dtype=int)
    else:
        filt_mask = np.load(filt_file)
    print(filt_mask)

    filters = load_filters(all_filters[filt_mask])
    #print(all_filters)
    #print(all_filters[filt_mask])
    obs['filters'] = filters
    obs['maggies'] = dat['phot_data'][gnum,z_idx][filt_mask]/3631
    obs['maggies_unc'] = dat['phot_data'][gnum,z_idx][filt_mask]/3631/sn
    obs['phot_mask'] = np.isfinite(dat['phot_data'][gnum,z_idx][filt_mask]/3631)
    obs['wavelength'] = None
    obs['spectrum'] = None
    obs['rest_sed'] = dat['sed_data'][gnum]/3631
    obs['wav'] = dat['wav_sed'][z_idx]
    obs['z_val'] = dat['z'][z_idx]
    obs['filt_config'] = filt_file
    return obs

def build_all(phot_file,galaxy_num,snr,z_idx,filt_file,**kwargs):

    return (build_obs(phot_file,galaxy_num,snr,z_idx,filt_file,**kwargs), build_model(phot_file,z_idx,**kwargs),
            build_sps(**kwargs))

run_params = {'verbose':False,
        'debug':False,
        'output_pickles': True,
        'nested_bound': 'multi', # bounding method
        'nested_sample': 'auto', # sampling method
        'nested_nlive_init': 400,
        'nested_nlive_batch': 200,
        'nested_bootstrap': 0,
        'nested_dlogz_init': 0.05,
        'nested_weight_kwargs': {"pfrac": 1.0},
        }


if __name__ == '__main__':

    # first version just run on all photometry
    # then have limited photometry setup
    # then factor in z term
    sim_id = str(sys.argv[1])
    snap = int(sys.argv[2])
    gal_num = int(sys.argv[3])
    full_phot = int(sys.argv[4])
    z_idx = int(sys.argv[5])
    filt_file = str(sys.argv[6])
    SNR = int(sys.argv[7])
    try:
        extra_label = str(sys.argv[8])
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

    outfile = f'{path_base}/prosp_runs/{sim_id}/snap{snap}/{prosp_out_folder}/prosp_run_{sim_id}_snap{snap}_galaxy{gal_num}{extra_label}.h5'
    print(outfile)
    #raise Exception()
    obs, model, sps = build_all(phot_file,gal_num,SNR,z_idx,filt_file,**run_params)
    
    run_params["sps_libraries"] = sps.ssp.libraries
    run_params["param_file"] = __file__
    print('Running fits')
    output = fit_model(obs, model, sps, [None,None],**run_params)
    print('Done. Writing now')
    print(model)
    print(obs)
    writer.write_hdf5(outfile, run_params, model, obs,
            output["sampling"][0], output["optimization"][0],
            tsample=output["sampling"][1],
            toptimize=output["optimization"][1])


import numpy as np
import yt
import caesar
import sys,os
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
from astropy import units as u
#gal = int(sys.argv[2])

def sed_load(filename):
        m = ModelOutput(filename)
        wav,flux = m.get_sed(inclination='all',aperture=-1)
        wav  = np.asarray(wav)#*u.micron #wav is in micron
        #wav *= (1.+z)
        flux = np.asarray(flux)[0]*u.erg/u.s
        return wav,flux

snap = sys.argv[1]

sim_id = "SIMBA_1P_2_0"

base_path = '/orange/narayanan/d.zimmerman/camels_test/'
caes_path = f'{base_path}catalogs_saved/{sim_id}/caesar_{snap}.hdf5'
pd_run_base_path = f'{base_path}pd_runs/{sim_id}/snap{snap}'

caes_file = caesar.load(caes_path)
valid_gals = np.loadtxt(f'{base_path}filtered/{sim_id}/snap{snap}/snap{snap}_gas_gals.txt',dtype=int)

gal_id = []
gmass_grid_l = []
gmass_grid_h = []
gmass_part = []
dmass_grid_l = []
dmass_grid_h = []
dmass_part = []


sed_low_old = []
sed_low_new = []
sed_high_old = []
sed_high_new = []

nref_low = 1
nref_high = 32

pd_run_low_old = f'{pd_run_base_path}_nref{nref_low}/'
pd_run_low_new = f'{pd_run_base_path}_nref{nref_low}_2/'

pd_run_high_old = f'{pd_run_base_path}_nref{nref_high}/'
pd_run_high_new = f'{pd_run_base_path}_nref{nref_high}_2/'


for gal in valid_gals: 
    print()
    print(gal)
    #wavl,sed_l = sed_load(f'{pd_run_low}snap{snap}.galaxy{gal}.rtout.sed')
    try:
        grid_prop_low = f'{pd_run_low_old}grid_physical_properties.{snap}_galaxy{gal}.npz'
        grid_prop_high = f'{pd_run_high_new}grid_physical_properties.{snap}_galaxy{gal}.npz'
        grid_data_l = np.load(grid_prop_low)
        grid_data_h = np.load(grid_prop_high)
        #print('h1')
        wavl,sed_l = sed_load(f'{pd_run_low_old}snap{snap}.galaxy{gal}.rtout.sed')
        wavh,sed_h = sed_load(f'{pd_run_high_old}snap{snap}.galaxy{gal}.rtout.sed')
        wavl2,sed_l2 = sed_load(f'{pd_run_low_new}snap{snap}.galaxy{gal}.rtout.sed')
        wavh2,sed_h2 = sed_load(f'{pd_run_high_new}snap{snap}.galaxy{gal}.rtout.sed')
        #print('h2')
    except:
        print(f'FAIL on {gal}')
        continue
    print()

    gal_id.append(gal)
    gmass_grid_l.append(np.sum(grid_data_l['grid_gas_mass']))
    gmass_grid_h.append(np.sum(grid_data_h['grid_gas_mass']))
    gmass_part.append(np.sum(grid_data_l['particle_gas_mass']))
    
    dmass_grid_l.append(np.sum(grid_data_l['grid_dustmass']))
    dmass_grid_h.append(np.sum(grid_data_h['grid_dustmass']))
    dmass_part.append(np.sum(grid_data_l['particle_dustmass']))


    sed_low_old.append(sed_l)
    sed_high_old.append(sed_h)

    sed_low_new.append(sed_l2)
    sed_high_new.append(sed_h2)

gal_id = np.array(gal_id)
gmass_grid_l = np.array(gmass_grid_l)
gmass_grid_h = np.array(gmass_grid_h)
gmass_part = np.array(gmass_part)

dmass_grid_l = np.array(dmass_grid_l)
dmass_grid_h = np.array(dmass_grid_h)
dmass_part = np.array(dmass_part)

sed_low_old = np.array(sed_low_old)
sed_high_old = np.array(sed_high_old)

sed_low_new = np.array(sed_low_new)
sed_high_new = np.array(sed_high_new)

'''
plt.scatter(gmass_part,gmass_grid_l/gmass_grid_h)
plt.xlabel('particle gas mass')
plt.ylabel('grid nref1/nref32 gas')
plt.title(f"z={np.round(caes_file.simulation.redshift,2)}")
plt.xscale('log')
plt.yscale('log')
plt.savefig(f'/home/d.zimmerman/figures/test_pd_gas_grid_{snap}_v4.png')
plt.close()

wav_repeat = np.repeat([wavl],len(gmass_grid_l),axis=0)

#for i in range(len(sed_low)):
#    plt.plot(wavl,sed_low[i]/sed_high[i],alpha=0.1)
plt.plot(wavl,np.median(sed_low/sed_high,axis=0),color = 'black')
plt.fill_between(wavl,np.quantile(sed_low/sed_high,0.16,axis=0),np.quantile(sed_low/sed_high,0.84,axis=0),alpha=0.5,color='grey')
plt.fill_between(wavl,np.quantile(sed_low/sed_high,0.05,axis=0),np.quantile(sed_low/sed_high,0.95,axis=0),alpha=0.2,color='grey')
plt.xlabel('wavelength')
plt.ylabel(f'sed nref{nref_low} (old) /nref{nref_high} (new) flux ratio')
plt.title(f"z={np.round(caes_file.simulation.redshift,2)}")
plt.ylim(10**-1,10**1)
plt.xscale('log')
plt.yscale('log')
plt.savefig(f'/home/d.zimmerman/figures/test_pd_sed_rat_{snap}_v4.png')
plt.close()
'''

#print(dmass_grid_h)
#print(dmass_part)

for i in range(len(gal_id)):
    gal = gal_id[i]
    if(np.abs(np.log10(dmass_grid_l[i])-np.log10(dmass_part[i])) < 1):
        low_plot = True
    else:
        low_plot = False

    if(np.abs(np.log10(dmass_grid_h[i])-np.log10(dmass_part[i])) < 1):
        high_plot = True
    else:
        high_plot = False

    if(not low_plot and not high_plot):
        continue
    print(gal)
    plt.plot(wavl,sed_low_new[i],label = 'new dust nref1',alpha=0.9)
    plt.plot(wavl,sed_high_new[i],label = 'new dust nref32',alpha=0.9)
    if(low_plot):
        plt.plot(wavl,sed_low_old[i],color = 'green',label = 'old dust nref1',alpha=0.7)
    if(high_plot):
        plt.plot(wavl,sed_high_old[i],color = 'magenta',label = 'old dust nref32',alpha=0.7)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('wavelength')
    plt.ylabel('flux')
    plt.legend()
    plt.title(f'galaxy {gal}')
    plt.savefig(f'/home/d.zimmerman/figures/pd_dust_comp/snap33_gal{gal}_close_dust_comp.png')
    plt.close()


import numpy as np
import matplotlib.pyplot as plt
import caesar
import yt 
snap = 33
sim_id = "SIMBA_1P_2_0"

base_path = '/orange/narayanan/d.zimmerman/camels_test/'
caes_path = f'{base_path}catalogs_saved/{sim_id}/caesar_{snap}.hdf5'
pd_run_base_path = f'{base_path}pd_runs/{sim_id}/snap{snap}'

caes_file = caesar.load(caes_path)
valid_gals = np.loadtxt(f'{base_path}filtered/{sim_id}/snap{snap}/snap{snap}_gas_gals.txt',dtype=int)

gmass_grid_l = []
gmass_grid_h = []
gmass_part = []
dmass_grid_l = []
dmass_grid_h = []
dmass_part = []
gal_id = []

nref_low = 1
nref_high = 32

pd_run_low = f'{pd_run_base_path}_nref{nref_low}/'
pd_run_high = f'{pd_run_base_path}_nref{nref_high}/'

for gal in valid_gals:
    print()
    print(gal)
    #wavl,sed_l = sed_load(f'{pd_run_low}snap{snap}.galaxy{gal}.rtout.sed')
    try:
        grid_prop_low = f'{pd_run_low}grid_physical_properties.{snap}_galaxy{gal}.npz'
        grid_prop_high = f'{pd_run_high}grid_physical_properties.{snap}_galaxy{gal}.npz'
        grid_data_l = np.load(grid_prop_low)
        grid_data_h = np.load(grid_prop_high)
    except:
        print(f'FAIL on {gal}')
        continue
    gal_id.append(gal)
    gmass_grid_l.append(np.sum(grid_data_l['grid_gas_mass']))
    gmass_grid_h.append(np.sum(grid_data_h['grid_gas_mass']))
    gmass_part.append(np.sum(grid_data_l['particle_gas_mass']))
    dmass_grid_l.append(np.sum(grid_data_l['grid_dustmass']))
    dmass_grid_h.append(np.sum(grid_data_h['grid_dustmass']))
    dmass_part.append(np.sum(grid_data_l['particle_dustmass']))

gal_id = np.array(gal_id)
gmass_grid_l = np.array(gmass_grid_l)
gmass_grid_h = np.array(gmass_grid_h)
gmass_part = np.array(gmass_part)
dmass_grid_l = np.array(dmass_grid_l)
dmass_grid_h = np.array(dmass_grid_h)
dmass_part = np.array(dmass_part)

arr = np.logspace(-3,3,1000)

plt.plot(arr,arr,color='black')
plt.scatter(gmass_grid_l/gmass_part,dmass_grid_l/dmass_part,label = f'{nref_low}',alpha=0.3)
plt.scatter(gmass_grid_h/gmass_part,dmass_grid_h/dmass_part,label = f'{nref_high}',alpha=0.3)

plt.xlabel('gas grid/particle gas mass')
plt.ylabel('dust grid/dust particle mass')
plt.title(f"z={np.round(caes_file.simulation.redshift,2)}")
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig(f'/home/d.zimmerman/figures/test_pd_gas_dust_grid_{snap}_v1.png')
plt.close()

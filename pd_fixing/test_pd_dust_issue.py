import numpy as np
import yt
import caesar
import sys,os
import matplotlib.pyplot as plt
#gal = int(sys.argv[2])
snap = sys.argv[1]

sim_id = "SIMBA_1P_2_0"

base_path = '/orange/narayanan/d.zimmerman/camels_test/'
caes_path = f'{base_path}catalogs_saved/{sim_id}/caesar_{snap}.hdf5'
pd_run_path = f'{base_path}pd_runs/{sim_id}/snap{snap}/'

caes_file = caesar.load(caes_path)
valid_gals = np.loadtxt(f'{base_path}filtered/{sim_id}/snap{snap}/snap{snap}_gas_gals.txt',dtype=int)

gal_id = []
ngas = []
smass = []
gmass = []
dmass_caes = []
dmass_yt = []
dmass_pd_part = []
dmass_pd_grid = []
ghm_rad = []
gmass_grid = []
gmass_part = []

for gal in valid_gals: 
#[560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578]:#valid_gals[358:508]:
    print()
    print(gal)
    try:
        filt_path = f'{base_path}filtered/{sim_id}/snap{snap}/galaxy_{gal}.hdf5'
        grid_prop = f'{pd_run_path}grid_physical_properties.{snap}_galaxy{gal}.npz'
        yt_data = yt.load(filt_path)
        filt_data = yt_data.all_data()
        grid_data = np.load(grid_prop)
    except:
        print(f'FAIL on {gal}')
        continue
    print()
    filt_dmass = np.log10(np.sum(filt_data[('PartType0','Dust_Masses')]))+np.log10(filt_data[('PartType0','Masses')][0].to('Msun').value/filt_data[('PartType0','Masses')][0].to('code_mass').value)
    part_dmass = np.log10(np.sum(grid_data['particle_dustmass']))
    grid_dmass = np.log10(np.sum(grid_data['grid_dustmass']))
    print('stellar:',np.log10(caes_file.galaxies[gal].masses['stellar'].to('Msun')))
    print('gas:',np.log10(caes_file.galaxies[gal].masses['gas'].to('Msun')))
    print('dust:',np.log10(caes_file.galaxies[gal].masses['dust'].to('Msun')))
    print('nstar:',caes_file.galaxies[gal].nstar)
    print('ngas:',caes_file.galaxies[gal].ngas)
    print()
    
    print('yt dust mass:',filt_dmass)
    print('pd particle dust mass:',part_dmass)
    print('pd grid dust mass:',grid_dmass)

    gal_id.append(gal)
    ngas.append(caes_file.galaxies[gal].ngas)
    smass.append(caes_file.galaxies[gal].masses['stellar'].to('Msun').value)
    gmass.append(caes_file.galaxies[gal].masses['gas'].to('Msun').value)
    dmass_caes.append(caes_file.galaxies[gal].masses['dust'].to('Msun').value)
    dmass_yt.append(np.sum(filt_data[('PartType0','Dust_Masses')])*filt_data[('PartType0','Masses')][0].to('Msun').value/filt_data[('PartType0','Masses')][0].to('code_mass').value)
    dmass_pd_part.append(np.sum(grid_data['particle_dustmass']))
    dmass_pd_grid.append(np.sum(grid_data['grid_dustmass']))
    ghm_rad.append(caes_file.galaxies[gal].radii['gas_r80'].to('kpc').value)
    gmass_grid.append(np.sum(grid_data['grid_gas_mass'])/2/10**30)
    gmass_part.append(np.sum(grid_data['particle_gas_mass'])/2/10**30)

gal_id = np.array(gal_id)
ngas = np.array(ngas)
smass = np.array(smass)
gmass = np.array(gmass)
dmass_caes = np.array(dmass_caes)
dmass_yt = np.array(dmass_yt)
dmass_pd_part = np.array(dmass_pd_part)
dmass_pd_grid = np.array(dmass_pd_grid)
ghm_rad = np.array(ghm_rad)
gmass_grid = np.array(gmass_grid)
gmass_part = np.array(gmass_part)

dmass_miss = dmass_yt < 1
dmass_grid = dmass_pd_grid < 1

oto = np.linspace(0,13,1000)

print(dmass_caes/gmass)
'''
plt.plot(oto,oto)
plt.scatter(np.log10(dmass_yt+1),np.log10(dmass_pd_grid+1),c=np.log10(smass))
plt.colorbar(label = 'log M*')
plt.xlabel('log+1 yt dust')
plt.ylabel('log+1 pd grid dust')
plt.title(f"z={np.round(caes_file.simulation.redshift,2)}")
plt.savefig(f'/home/d.zimmerman/figures/test_dustgrid_{snap}_v1.png')
plt.close()

#plt.scatter(np.log10(smass_yt),np.log10(dmass_yt+1),c=np.log10(dmass_yt/gmass+1),label='good')
plt.scatter(np.log10(smass),np.log10(dmass_yt+1),label='good')
plt.scatter(np.log10(smass)[dmass_grid],np.log10(dmass_yt+1)[dmass_grid],label='grid bad')
plt.scatter(np.log10(smass)[dmass_miss],np.log10(dmass_yt+1)[dmass_miss],label='gas but no dust')
plt.xlabel('log M*')
plt.ylabel('yt log Md+1')
plt.legend()
plt.title(f"z={np.round(caes_file.simulation.redshift,2)}")
plt.savefig(f'/home/d.zimmerman/figures/test_dustgrid_cat_check_{snap}_vdmass.png')
plt.close()


plt.scatter(np.log10(smass),ghm_rad,label='good')
plt.scatter(np.log10(smass)[dmass_grid],ghm_rad[dmass_grid],label='grid bad')
plt.scatter(np.log10(smass)[dmass_miss],ghm_rad[dmass_miss],label='gas but no dust')
plt.xlabel('log M*')
plt.ylabel('gas 80% radius')
plt.legend()
plt.title(f"z={np.round(caes_file.simulation.redshift,2)}")
plt.savefig(f'/home/d.zimmerman/figures/test_dustgrid_cat_check_{snap}_vhmr.png')
plt.close()



plt.scatter(np.log10(smass),np.log10(dmass_yt/gmass+1),label='good')
plt.scatter(np.log10(smass)[dmass_grid],np.log10(dmass_yt/gmass+1)[dmass_grid],label='grid bad')
plt.scatter(np.log10(smass)[dmass_miss],np.log10(dmass_yt/gmass+1)[dmass_miss],label='gas but no dust')
plt.xlabel('log M*')
plt.ylabel('log dtg+1')
plt.legend()
plt.title(f"z={np.round(caes_file.simulation.redshift,2)}")
plt.savefig(f'/home/d.zimmerman/figures/test_dustgrid_cat_check_{snap}_vdtg.png')
plt.close()

plt.scatter(np.log10(smass),np.log10((dmass_yt+1)/gmass),label='good')
plt.scatter(np.log10(smass)[dmass_grid],np.log10((dmass_yt+1)/gmass)[dmass_grid],label='grid bad')
plt.scatter(np.log10(smass)[dmass_miss],np.log10((dmass_yt+1)/gmass)[dmass_miss],label='gas but no dust')
plt.xlabel('log M*')
plt.ylabel('log dtg+1')
plt.legend()
plt.title(f"z={np.round(caes_file.simulation.redshift,2)}")
plt.savefig(f'/home/d.zimmerman/figures/test_dustgrid_cat_check_{snap}_vdtg.png')
plt.close()

plt.scatter(np.log10(gmass/smass),ghm_rad,label='good')
plt.scatter(np.log10(gmass/smass)[dmass_grid],ghm_rad[dmass_grid],label='grid bad')
plt.scatter(np.log10(gmass/smass)[dmass_miss],ghm_rad[dmass_miss],label='gas but no dust')
plt.xlabel('log gts')
plt.ylabel('gas 80% rad')
plt.legend()
plt.title(f"z={np.round(caes_file.simulation.redshift,2)}")
plt.savefig(f'/home/d.zimmerman/figures/test_dustgrid_cat_check_{snap}_vgts.png')
plt.close()
'''

plt.plot(oto,oto)
plt.scatter(np.log10((gmass_part+1)),np.log10((gmass_grid+1)),c=np.log10(smass))
plt.colorbar(label = 'log M*')
plt.xlabel('log+1 particle gas')
plt.ylabel('log+1 grid gas')
plt.title(f"z={np.round(caes_file.simulation.redshift,2)}")
plt.savefig(f'/home/d.zimmerman/figures/test_gasgrid_{snap}_v1.png')
plt.close()

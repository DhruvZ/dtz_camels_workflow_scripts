import yt
import numpy as np

use_gather = True

nref_array = [1024,512,256,128,64,32,16,8]
ds = yt.load('gizmo_mhd_mwdisk/gizmo_mhd_mwdisk.hdf5')
sm = 'scat'
if(use_gather):
    ds.sph_smoothing_style = "gather"
    sm = 'gath'

ad = ds.all_data()
part_gas_mass = np.sum(ad[('PartType0','Masses')].to('Msun').value)


n_ratios = []
n_masses = []
mass_bad = []
for nval in nref_array:
    octree = ds.octree(ds.domain_left_edge,ds.domain_right_edge,n_ref=nval)
    valb = np.sum(octree["PartType0","Masses"].to('Msun').value)
    val0 = octree["PartType0","Density"].to('Msun/kpc**3')
    #q = octree["PartType0","Masses"].to('Msun')
    refined = octree[('index','refined')].astype('bool')
    volume = octree['index','dx'][~refined]*octree['index','dy'][~refined]*octree['index','dz'][~refined]
    val = np.sum(volume.to('kpc**3')*val0)
    #print(octree['index','dx'][~refined])
    #print(len(octree['index','dx'][~refined]))
    #print(len(octree["PartType0","Density"].to('Msun/kpc**3').value))
    mass_bad.append(valb)
    n_masses.append(val)
    n_ratios.append(val/part_gas_mass)

n_masses = np.array(n_masses)
n_ratios = np.array(n_ratios)
mass_bad = np.array(mass_bad)
print('particle')
print(part_gas_mass)
print('octree')
print(nref_array)
print(n_masses)
print(n_ratios)
print('bad')
print(mass_bad)
print(part_gas_mass/mass_bad)
np.save(f'mass_ratios_{sm}.npy',n_ratios)

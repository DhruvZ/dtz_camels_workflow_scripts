import numpy as np
import yt
import h5py
import yt.units as u
import matplotlib.pyplot as plt

snap_num = 33
good_gals = np.loadtxt(f'/orange/narayanan/d.zimmerman/camels_test/filtered/SIMBA_1P_2_0/snap{snap_num}/snap{snap_num}_gas_gals.txt',dtype=int)
nref_checks = [1,2,4,8,16,32]
bbox_lim = 1.e5

gal_pos = np.load(f'/orange/narayanan/d.zimmerman/camels_test/filtered/SIMBA_1P_2_0/snap{snap_num}/snap{snap_num}_gal_positions.npz',allow_pickle=True)
print(gal_pos['ngalaxies'])
all_pos = gal_pos['pos'][()]

box_len0 = 100

oct_m = []
part_m = []


def get_masses(fname,pos):
    x_cent = pos[0]#19496.465
    y_cent = pos[1]#9300.578
    z_cent = pos[2]#6220.7397

    ds = yt.load(fname)#,bounding_box = bbox)

    center = [x_cent,y_cent,z_cent]
    print ('[octree zoom_bbox_filter:] using center: ',center)
    box_len = box_len0
    box_len = ds.quan(box_len,'kpc')
    box_len = float(box_len.to('code_length').value)
    bbox_lim = box_len
    bbox2 = [[center[0]-bbox_lim,center[0]+bbox_lim],
            [center[1]-bbox_lim,center[1]+bbox_lim],
            [center[2]-bbox_lim,center[2]+bbox_lim]]
    print()
    print ('[octree zoom] new zoomed bbox (comoving/h) in code units= ',bbox2)

    gal_octs = []

    left = np.array([pos[0] for pos in bbox2])
    right = np.array([pos[1] for pos in bbox2])
    reg = ds.all_data()
    part_mass = np.sum(reg[("PartType0","Masses")].to('Msun'))
    for nval in nref_checks:
        octree = ds.octree(left,right,n_ref=nval)
        #oct_mass = np.sum(octree["PartType0","Masses"].to('Msun'))
        val0 = octree["PartType0","Density"].to('Msun/kpc**3')
        refined = octree[('index','refined')].astype('bool')
        volume = octree['index','dx'][~refined]*octree['index','dy'][~refined]*octree['index','dz'][~refined]
        oct_mass = np.sum(volume.to('kpc**3')*val0)
        gal_octs.append(oct_mass)
    ds.close()
    oct_m.append(gal_octs)
    
    
    return part_mass

for galaxy_num in good_gals:
    print()
    print()
    print('gal:',galaxy_num)
    print()
    fname = f'/orange/narayanan/d.zimmerman/camels_test/filtered/SIMBA_1P_2_0/snap{snap_num}/galaxy_{galaxy_num}.hdf5'# name of filtered file
    pos = all_pos[f'galaxy{galaxy_num}'][f'snap{snap_num}']
    print(pos)
    pm = get_masses(fname,pos)
    part_m.append(pm)



oct_m = np.array(oct_m)
part_m = np.array(part_m)

for i in range(len(nref_checks)):
    plt.scatter(part_m,oct_m[:,i]/part_m,label=f'nref{nref_checks[i]}',alpha=0.3)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Mg (particle)')
plt.ylabel('Mg smoothed/Mg particle')
plt.legend()
plt.savefig(f'/home/d.zimmerman/figures/nref_eval_snap{snap_num}.png')



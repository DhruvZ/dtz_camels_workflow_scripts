import numpy as np
import yt
import h5py
import yt.units as u
import matplotlib.pyplot as plt

snap_num = 33
good_gals = np.loadtxt(f'/orange/narayanan/d.zimmerman/camels_test/filtered/SIMBA_1P_2_0/snap{snap_num}/snap{snap_num}_gas_gals.txt',dtype=int)
nref_low = 1
nref_high = 16
bbox_lim = 1.e5

gal_pos = np.load(f'/orange/narayanan/d.zimmerman/camels_test/filtered/SIMBA_1P_2_0/snap{snap_num}/snap{snap_num}_gal_positions.npz',allow_pickle=True)
print(gal_pos['ngalaxies'])
all_pos = gal_pos['pos'][()]

box_len0 = 100
#    bbox2 = [[-2.*bbox_lim,2.*bbox_lim],
#            [-2.*bbox_lim,2.*bbox_lim],
#            [-2.*bbox_lim,2.*bbox_lim]]

oct_low_m = []
oct_high_m = []
part_low_m = []
part_high_m = []


def get_masses(fname,nref,pos):
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

    left = np.array([pos[0] for pos in bbox2])
    right = np.array([pos[1] for pos in bbox2])
    octree = ds.octree(left,right,n_ref=nref)
    reg = ds.all_data()
    print(octree)
    oct_mass = np.sum(octree["PartType0","Masses"].to('Msun'))
    part_mass = np.sum(reg[("PartType0","Masses")].to('Msun'))
    ds.close()
    return oct_mass,part_mass

for galaxy_num in good_gals:
    print()
    print()
    print('gal:',galaxy_num)
    print()
    fname = f'/orange/narayanan/d.zimmerman/camels_test/filtered/SIMBA_1P_2_0/snap{snap_num}/galaxy_{galaxy_num}.hdf5'# name of filtered file
    pos = all_pos[f'galaxy{galaxy_num}'][f'snap{snap_num}']
    print(pos)
    low_oct, low_part = get_masses(fname,nref_low,pos)
    high_oct, high_part = get_masses(fname,nref_high,pos)

    oct_low_m.append(low_oct)
    oct_high_m.append(high_oct)
    part_low_m.append(low_part)
    part_high_m.append(high_part)



oct_low_m = np.array(oct_low_m)
oct_high_m = np.array(oct_high_m)
part_m = np.array(part_low_m)

plt.scatter(part_m,oct_low_m/oct_high_m)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Mg (particle)')
plt.ylabel('Mg nref = 1/nref = 16')
plt.savefig(f'/home/d.zimmerman/figures/nref_dependence_snap{snap_num}.png')



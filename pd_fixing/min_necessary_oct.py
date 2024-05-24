import numpy as np
import yt
import h5py
import yt.units as u
import pdb
galaxy_num = 78
fname = f'/orange/narayanan/d.zimmerman/camels_test/filtered/SIMBA_1P_2_0/snap33/galaxy_{galaxy_num}.hdf5'# name of filtered file
nref = 1
bbox_lim = 1.e5
gal_pos = np.load(f'/orange/narayanan/d.zimmerman/camels_test/filtered/SIMBA_1P_2_0/snap33/snap33_gal_positions.npz',allow_pickle=True)
print(gal_pos['ngalaxies'])
pos = gal_pos['pos'][()]
pos = pos[f'galaxy{galaxy_num}']['snap33']
print(pos)
x_cent = pos[0]#19496.465
y_cent = pos[1]#9300.578
z_cent = pos[2]#6220.7397

center = [x_cent,y_cent,z_cent]
print ('[octree zoom_bbox_filter:] using center: ',center)

bbox_max = False
gather = True

ds = yt.load(fname)#,bounding_box = bbox)
#ds.index
if(bbox_max):
    bbox = [[-2.*bbox_lim,2.*bbox_lim],
            [-2.*bbox_lim,2.*bbox_lim],
            [-2.*bbox_lim,2.*bbox_lim]]
else:
    box_len = 100
    box_len = ds.quan(box_len,'kpc')
    box_len = float(box_len.to('code_length').value)
    bbox_lim = box_len
    bbox = [[center[0]-bbox_lim,center[0]+bbox_lim],
            [center[1]-bbox_lim,center[1]+bbox_lim],
            [center[2]-bbox_lim,center[2]+bbox_lim]]
    print ('[octree zoom] new zoomed bbox (comoving/h) in code units= ',bbox)


if(gather):
    ds.sph_smoothing_style = "gather"
print()

left = np.array([pos[0] for pos in bbox])
right = np.array([pos[1] for pos in bbox])
#breakpoint()
octree = ds.octree(left,right,n_ref=nref)
reg = ds.all_data()
mass = octree["PartType0","Masses"].to('Msun')
print(octree)
#breakpoint()
print(octree["PartType0","Masses"].to('Msun')/10**6)
#print(reg[("gas","smoothedmasses")].to('Msun'))
print(reg[("gas","mass")].to('Msun')/10**6)
print(np.log10(np.sum(octree["PartType0","Masses"].to('Msun'))))
print(np.log10(np.sum(reg[("PartType0","Masses")].to('Msun'))))

import numpy as np
import yt
import h5py
import yt.units as u
import matplotlib.pyplot as plt

galaxy_num = 560
fname = f'/orange/narayanan/d.zimmerman/camels_test/filtered/SIMBA_1P_2_0/snap33/galaxy_{galaxy_num}.hdf5'# name of filtered file
nref = 16
bbox_lim = 1.e5
gal_pos = np.load(f'/orange/narayanan/d.zimmerman/camels_test/filtered/SIMBA_1P_2_0/snap33/snap33_gal_positions.npz',allow_pickle=True)
print(gal_pos['ngalaxies'])
pos = gal_pos['pos'][()]
pos = pos[f'galaxy{galaxy_num}']['snap33']
print(pos)
x_cent = pos[0]#19496.465
y_cent = pos[1]#9300.578
z_cent = pos[2]#6220.7397

bbox2 = [[-2.*bbox_lim,2.*bbox_lim],
         [-2.*bbox_lim,2.*bbox_lim],
         [-2.*bbox_lim,2.*bbox_lim]]
ds = yt.load(fname)#,bounding_box = bbox)

center = [x_cent,y_cent,z_cent]
print ('[octree zoom_bbox_filter:] using center: ',center)
box_len = 50
box_len = ds.quan(box_len,'kpc')
box_len = float(box_len.to('code_length').value)
bbox_lim = box_len
bbox2 = [[center[0]-bbox_lim,center[0]+bbox_lim],
            [center[1]-bbox_lim,center[1]+bbox_lim],
            [center[2]-bbox_lim,center[2]+bbox_lim]]
#center = ds.quan(center,'code_length')
print()
print ('[octree zoom] new zoomed bbox (comoving/h) in code units= ',bbox2)

left = np.array([pos[0] for pos in bbox2])
right = np.array([pos[1] for pos in bbox2])
octree = ds.octree(left,right,n_ref=nref)
reg = ds.all_data()

print(octree)
print(octree["PartType0","Masses"].to('Msun'))
print(reg[("PartType0","Masses")].to('Msun'))
print(np.log10(np.sum(octree["PartType0","Masses"].to('Msun'))))
print(np.log10(np.sum(reg[("PartType0","Masses")].to('Msun'))))



prj = yt.ParticlePlot(ds,('PartType0','particle_position_x'),('PartType0','particle_position_y'),
        ('PartType0','mass'),center=center,width=(box_len,box_len))
prj.save(f"/home/d.zimmerman/figures/test_z0_gal_{galaxy_num}_filt")

prj = yt.ParticlePlot(ds,('PartType0','particle_position_x'),('PartType0','particle_position_y'),
        ('PartType0','density'),center=center,width=(box_len,box_len))
prj.save(f"/home/d.zimmerman/figures/test_z0_gal_{galaxy_num}_filt")


prj = yt.ProjectionPlot(ds,"z",("PartType0","density"),center=pos,width = box_len)
prj.save(f"/home/d.zimmerman/figures/test_z0_gal_{galaxy_num}_filt")



#prj2 = yt.ParticlePlot(yt_filt,('all','particle_position_x'),('all','particle_position_y'),('all','particle_mass'),width=(wid,wid))

import h5py
import numpy as np
import sys, os
import tqdm
import pandas as pd
import yt
###############
# Line arguments
##############
sim = sys.argv[1]
run_1p_param = int(sys.argv[2])
run_1p_num = sys.argv[3]
snap_num = int(sys.argv[4])
max_num = int(sys.argv[5])

snapshot_path = f'/orange/narayanan/d.zimmerman/camels_results/sims_loaded/{sim}_1P_p'

sim_id = f'{sim}_1P_p{run_1p_param}_{run_1p_num}'
output_path = f'/orange/narayanan/d.zimmerman/camels_results/filtered/{sim_id}/snap{snap_num}/'
ds = f"{snapshot_path}{run_1p_param}_{run_1p_num}/snapshot_"

if(int(snap_num) < 10):
    snap_str = f"00{snap_num}"
elif(int(snap_num) < 100):
    snap_str = f"0{snap_num}"
else:
    snap_str = f"{snap_num}"

input_file = h5py.File(f'{ds}{snap_str}.hdf5', 'r')

cat_file = f'/orange/narayanan/d.zimmerman/camels_results/catalogs_saved/{sim_id}/groups_{snap_num}.hdf5'
catalog = h5py.File(cat_file,'r')


group_len_gas = catalog['Group/GroupLenType'][:,0]
group_len_star = catalog['Group/GroupLenType'][:,4]

subhalo_len_gas = catalog['Subhalo/SubhaloLenType'][:,0]
subhalo_len_star = catalog['Subhalo/SubhaloLenType'][:,4]

group_substart = catalog['Group/GroupFirstSub'][:]
group_subcount = catalog['Group/GroupNsubs'][:]
subhalo_groupnum = catalog['Subhalo/SubhaloGrNr'][:]


num = np.array(range(len(subhalo_len_gas)),dtype=int)

gal_list = num[subhalo_len_star >= 24]
gal_count = len(gal_list)

if(max_num < 0):
    gal_max = gal_count
else:
    gal_max = max_num


for galaxy in range(gal_max):
    subhalo = gal_list[int(galaxy)]
    outfile = f'{output_path}/galaxy_{galaxy}.hdf5'
    print('creating snapshot for galaxy '+str(galaxy))
    print('creating snapshot for subhalo '+str(subhalo))
    
    with h5py.File(outfile, 'w') as output_file: 
        output_file.copy(input_file['Header'], 'Header')
        output_file.copy(catalog['Config'], 'Config')
        
        print('copying gas attributes')
        output_file.create_group('PartType0')
        #print('/Offsets/'+str(snap_num)+'/Subhalo/SnapByType')
        #print('/Groups/'+str(snap_num)+'/Subhalo/SubhaloLenType')
        
        start = np.sum(group_len_gas[:subhalo_groupnum[subhalo]])+np.sum(subhalo_len_gas[group_substart[subhalo_groupnum[subhalo]]:subhalo])
        glength = subhalo_len_gas[subhalo]
        #print('gas start: ',start)
        #print('gas length: ',glength)
        #print('gas end: ',start+glength)
        for k in tqdm.tqdm(input_file['/PartType0/']):
            output_file['PartType0'][str(k)] = input_file['/PartType0/'+str(k)][start:start+glength]
        
        print('copying star attributes')
        output_file.create_group('PartType4')
        
        start = np.sum(group_len_star[:subhalo_groupnum[subhalo]])+np.sum(subhalo_len_star[group_substart[subhalo_groupnum[subhalo]]:subhalo])
        slength = subhalo_len_star[subhalo]
        
        #print('star start: ',start)
        #print('star length: ',slength)
        #print('star end: ',start+slength)
        for k in tqdm.tqdm(input_file['/PartType4/']):
            c = str(k)
            if c == 'StellarAssembly':
                continue
            output_file['PartType4'][str(k)] = input_file['/PartType4/'+str(k)][start:start+slength]

        
        
    re_out = h5py.File(outfile, 'r+')
        
    re_out['Header'].attrs.modify('NumFilesPerSnapshot', 1)
    re_out['Header'].attrs.modify('NumPart_ThisFile', np.array([glength, 0, 0, 0, slength, 0]))
    re_out['Header'].attrs.modify('NumPart_Total', np.array([glength, 0, 0, 0, slength, 0]))
        
    re_out.close()
        
    print('galaxy '+str(galaxy)+' is done.')
    print()

np.save(f'{output_path}/snap{snap_num}_subhaloIDs.npy',gal_list)


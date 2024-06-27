import numpy as np
import sys,os,yt
import caesar
import argparse


sim = str(sys.argv[1])
param_1p = int(sys.argv[2])
param_run = str(sys.argv[3])
snap = int(sys.argv[4])


if(snap < 10):
    snap_str = f'00{snap}'
elif(snap < 100):
    snap_str = f'0{snap}'
else:
    snap_str = f'{snap}'


sim_id = f'{sim}_1P_p{param_1p}_{param_run}'

base_path = '/orange/narayanan/d.zimmerman/camels_results'
sim_loc = f'{base_path}/sims_loaded/{sim_id}/snapshot_{snap_str}.hdf5'

###################################### INPUT ##########################################
#root           = '/mnt/ceph/users/camels/PUBLIC_RELEASE'
ssp_table_file = '/orange/narayanan/d.zimmerman/camels_scripts/FSPS_Chab_EL.hdf5'
nproc          = 8
#######################################################################################


f_fof    = f'{base_path}/catalogs_loaded/{sim_id}/fof6d2_{snap_str}'
fout     = f'{base_path}/catalogs_loaded/{sim_id}/caesar_newsnaps2_{snap_str}.hdf5'#%(folder_out,snapnum)

ds = yt.load(sim_loc)
        
obj = caesar.CAESAR(ds)
    
obj.member_search(haloid='fof', fof6d_file=f_fof, fsps_bands='uvoir', 
        ssp_model='FSPS', ssp_table_file=ssp_table_file, 
        ext_law='composite', nproc=nproc)

# save the CAESAR galaxy/halo catalogue
obj.save(fout)

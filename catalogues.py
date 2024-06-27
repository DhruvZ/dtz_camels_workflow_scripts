import numpy as np
import sys,os,yt
import caesar
import argparse

# read the elements to be removed and the ones already discarded
#parser = argparse.ArgumentParser(description="description of routine")
#parser.add_argument("-folder", default=[], help="simulation name")
#parser.add_argument("-sim", default=[], help="simulation name")
#args   = parser.parse_args()
#folder = args.folder
#sim    = args.sim


###################################### INPUT ##########################################
#root           = '/mnt/ceph/users/camels/PUBLIC_RELEASE'
ssp_table_file = 'FSPS_Chab_EL.hdf5'
nproc          = 10
#######################################################################################

# create output folder if it doesnt exist
folder_out = '%s/Caesar/%s/%s'%(root,sim,folder)
if not(os.path.exists(folder_out)):  os.system('mkdir %s'%folder_out)

# do a loop over all redshifts
for snapnum in range(90):

    # get the name of the files
    snapshot = '%s/Sims/%s/%s/snap_%03d.hdf5'%(root,sim,folder,snapnum)
    f_fof    = '%s/fof6d_%03d'%(folder_out,snapnum)
    fout     = '%s/caesar_%03d.hdf5'%(folder_out,snapnum)
    if os.path.exists(fout):  continue

    # load the snapshot into yt
    ds = yt.load(snapshot)
        
    # create a CAESAR object, and pass along the yt dataset
    obj = caesar.CAESAR(ds)
    
    # execute member_search(), which identifies halos, galaxies, and 
    # computes properties, including photometry
    obj.member_search(haloid='fof', fof6d_file=f_fof, fsps_bands='uvoir', 
                      ssp_model='FSPS', ssp_table_file=ssp_table_file, 
                      ext_law='composite', nproc=nproc)
    
    # save the CAESAR galaxy/halo catalogue
    obj.save(fout)

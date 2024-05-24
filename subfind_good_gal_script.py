import yt
import caesar
import numpy as np
import sys
import matplotlib.pyplot as plt
import h5py

sim_id = sys.argv[1]
num = int(sys.argv[2])

fileroot = f'/orange/narayanan/d.zimmerman/camels_results/catalogs_saved/{sim_id}/group_'
saveroot = f'/orange/narayanan/d.zimmerman/camels_results/filtered/{sim_id}/snap'
fileex='.hdf5'

gal_ids = np.load(f'{saveroot}{num}/snap{num}_subhaloIDs.npy',dtype=int)


#for num in snapnums:
group_obj = h5py.File(fileroot+str(num)+fileex,'r')

gal_gasses = group_obj['Subhalo/SubhaloMassType'][:,0][gal_ids]*10.0**10
gal_index_list = np.array(range(len(gal_gasses)),dtype=int)
#print(gal_index_list)
good_gals = gal_index_list[gal_gasses > 0]
#print(good_gals)
test_file = open(saveroot+str(num)+"/snap"+str(num)+"_gas_gals.txt","w")
#for j in good_gals:
np.savetxt(test_file,good_gals,fmt='%s')
test_file.close()

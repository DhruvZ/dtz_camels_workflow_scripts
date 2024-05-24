import yt
import caesar
import numpy as np
import sys
import matplotlib.pyplot as plt

sim_id = sys.argv[1]
num = int(sys.argv[2])

fileroot = f'/orange/narayanan/d.zimmerman/camels_results/catalogs_saved/{sim_id}/caesar_'
saveroot = f'/orange/narayanan/d.zimmerman/camels_results/filtered/{sim_id}/snap'
fileex='.hdf5'


#for num in snapnums:
caes_obj = caesar.load(fileroot+str(num)+fileex)
gal_gasses = np.array([caes_obj.galaxies[i].masses['gas'] for i in range(len(caes_obj.galaxies))])
gal_index_list = np.array(range(len(gal_gasses)),dtype=int)
#print(gal_index_list)
good_gals = gal_index_list[gal_gasses > 0]
#print(good_gals)
test_file = open(saveroot+str(num)+"/snap"+str(num)+"_gas_gals.txt","w")
#for j in good_gals:
np.savetxt(test_file,good_gals,fmt='%s')
test_file.close()

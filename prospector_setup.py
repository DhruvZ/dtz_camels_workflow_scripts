#purpose: to set up slurm files for prospector runs 
#for the DTZ CAMELS work.This is written for 
#the University of Florida's HiPerGator2 cluster.
import yt
import caesar
import numpy as np
from subprocess import call
import sys

#################
sim_id = str(sys.argv[1])
snap_num = int(sys.argv[2])

try:
    base_path = str(sys.argv[3])
except:
    base_path = '/orange/narayanan/d.zimmerman/camels_results/'


print(f'{base_path}/prosp_runs/{sim_id}/')

npzfile = f'{base_path}/ml_data/all_filter_{sim_id}_snap{snap_num}.npz'  
data = np.load(npzfile)
ngals = len(data['gal_num'])


#################

cmd = "/orange/narayanan/d.zimmerman/camels_scripts/camels_setup_all_prosp.hpg.sh "+str(base_path)+' '+str(sim_id)+' '+str(snap_num)+' '+str(ngals)
call(cmd,shell=True)

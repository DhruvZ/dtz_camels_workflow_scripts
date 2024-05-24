#purpose: to set up slurm files and model *.py files from the
#positions written by caesar_cosmology_npzgen.py for a cosmological
#simulation.  This is written for the University of Florida's
#HiPerGator2 cluster.
import yt
import caesar
import numpy as np
from subprocess import call
import sys

nnodes=1


#################
# Edit these !!!
sim_id = sys.argv[1]
snap_num = sys.argv[2]

caes_file = caesar.load(f'/orange/narayanan/d.zimmerman/camels_results/catalogs_saved/{sim_id}/caesar_{snap_num}.hdf5')
snap_redshift = caes_file.simulation.redshift
gal_count = len(caes_file.galaxies)

npzfile = f'/orange/narayanan/d.zimmerman/camels_results/filtered/{sim_id}/snap{snap_num}/snap{snap_num}_gal_positions.npz'  
model_dir_base = f'/orange/narayanan/d.zimmerman/camels_results/pd_scripts/{sim_id}/snap{snap_num}/'
out_dir_base = f'/orange/narayanan/d.zimmerman/camels_results/pd_runs/{sim_id}/snap{snap_num}/'
hydro_dir = f'/orange/narayanan/d.zimmerman/camels_results/filtered/{sim_id}/snap{snap_num}/'
hydro_dir_remote = hydro_dir

#model_run_name='simba_m100n1024'
model_run_name=f'{sim_id}'

#################

COSMOFLAG=0 #flag for setting if the gadget snapshots are broken up into multiples or not and follow a nomenclature snapshot_000.0.hdf5
FILTERFLAG = 1 #flag for setting if the gadget snapshots are filtered or not, and follow a nomenclature snap305_galaxy1800_filtered.hdf5


SPHGR_COORDINATE_REWRITE = True


#===============================================

if (COSMOFLAG == 1) and (FILTERFLAG == 1):
    raise ValueError("COSMOFLAG AND FILTER FLAG CAN'T BOTH BE SET")


data = np.load(npzfile,allow_pickle=True)
pos = data['pos'][()] #positions dictionary
#ngalaxies is the dict that says how many galaxies each snapshot has, in case it's less than NGALAXIES_MAX
ngalaxies = data['ngalaxies'][()]
print(ngalaxies)

#print("STARTING LOOP")
g_gal_count = len(np.loadtxt(f'/orange/narayanan/d.zimmerman/camels_results/filtered/{sim_id}/snap{snap_num}/snap{snap_num}_gas_gals.txt'))
for snap in [snap_num]:
    
    model_dir = model_dir_base#+'/snap{:03d}'.format(snap)
    model_dir_remote = model_dir
    
    tcmb = 2.73*(1.+snap_redshift)

    NGALAXIES = ngalaxies['snap'+str(snap)]
    #print("INITIAL DONE")
	
    for nh in range(NGALAXIES):
        try:
            xpos = pos['galaxy'+str(nh)]['snap'+str(snap)][0]
        except: continue
        
        ypos = pos['galaxy'+str(nh)]['snap'+str(snap)][1]
        zpos = pos['galaxy'+str(nh)]['snap'+str(snap)][2]
        #print("CALLING")	
        cmd = "/orange/narayanan/d.zimmerman/camels_scripts/camels_setup_all_pd.hpg.sh "+str(nnodes)+' '+model_dir+' '+hydro_dir+' '+out_dir_base+' '+model_run_name+' '+str(COSMOFLAG)+' '+str(FILTERFLAG)+' '+model_dir_remote+' '+hydro_dir_remote+' '+str(xpos)+' '+str(ypos)+' '+str(zpos)+' '+str(nh)+' '+str(snap)+' '+str(tcmb)+' '+str(g_gal_count-1)
        call(cmd,shell=True)
        #print("CALLED")
        

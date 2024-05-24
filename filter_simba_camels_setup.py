import h5py
import caesar
import sys
import glob
import numpy as np
import tqdm


###########
# Line arguments
###########
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

caesar_file = f'/orange/narayanan/d.zimmerman/camels_results/catalogs_saved/{sim_id}/caesar_{snap_num}.hdf5'
obj = caesar.load(caesar_file)

gal_count = len(obj.galaxies)

if(max_num < 0):
    gal_max = gal_count
else:
    gal_max = max_num

for galaxy in range(gal_max):
    print()
    print("GALAXY NUM:",str(galaxy))
    glist = obj.galaxies[int(galaxy)].glist
    slist = obj.galaxies[int(galaxy)].slist
    print('star count:',len(slist))
    print('gas count:',len(glist))
    print()
    
    with h5py.File(output_path+'galaxy_'+str(galaxy)+'.hdf5', 'w') as output_file:
        output_file.copy(input_file['Header'], 'Header')
        print('starting with gas attributes now')
        output_file.create_group('PartType0')
        for k in tqdm.tqdm(input_file['PartType0']):
            output_file['PartType0'][k] = input_file['PartType0'][k][:][glist]
        print('moving to star attributes now')
        output_file.create_group('PartType4')
        for k in tqdm.tqdm(input_file['PartType4']):
            output_file['PartType4'][k] = input_file['PartType4'][k][:][slist]
        
        print('done copying attributes, going to edit header now')
        outfile_reload = output_path+'galaxy_'+str(galaxy)+'.hdf5'
        
        re_out = h5py.File(outfile_reload,'r+')                                                  
        re_out['Header'].attrs.modify('NumPart_ThisFile', np.array([len(glist), 0, 0, 0, len(slist), 0]))  
        re_out['Header'].attrs.modify('NumPart_Total', np.array([len(glist), 0, 0, 0, len(slist), 0]))  
        
        re_out.close()

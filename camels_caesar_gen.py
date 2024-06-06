import yt
import caesar
import sys,os

sim = str(sys.argv[1])
run_1p_param = int(sys.argv[2])
run_1p_num = str(sys.argv[3])
snap_num = int(sys.argv[4])

sim_id = f'{sim}_1P_p{run_1p_param}_{run_1p_num}'
sim_dir = f'/orange/narayanan/d.zimmerman/camels_results/sims_loaded/{sim_id}/'
pro = 16

if(snap_num < 10):
    snap_str = f'00{snap_num}'
elif(snap_num < 100):
    snap_str = f'0{snap_num}'
else:
    snap_str = f'{snap_num}'
    

caesar_save_loc = f'/orange/narayanan/d.zimmerman/camels_results/catalogs_loaded/{sim_id}/caesar_newsnaps_{snap_str}.hdf5'

ds = yt.load(f'{sim_dir}/snapshot_{snap_str}.hdf5')
print()
print("loaded")
print()
caes_obj = caesar.CAESAR(ds)
print()
print("object setup")
print()
caes_obj.member_search(haloid = 'fof',blackholes=True,skipran=False,fof6d = True,nproc = pro)
print()
print("member search")
print()
caes_obj.save(caesar_save_loc)

print()
print("saved")
print()



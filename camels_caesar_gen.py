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

ssp_table_file = 'FSPS_Chab_EL.hdf5'

if(snap_num < 10):
    snap_str = f'00{snap_num}'
elif(snap_num < 100):
    snap_str = f'0{snap_num}'
else:
    snap_str = f'{snap_num}'
    

caesar_save_loc = f'/orange/narayanan/d.zimmerman/camels_results/catalogs_loaded/{sim_id}/caesar_newsnaps_{snap_str}.hdf5'
fof_save_loc = f'/orange/narayanan/d.zimmerman/camels_results/catalogs_loaded/{sim_id}/fof6d_newsnaps_{snap_str}'

ds = yt.load(f'{sim_dir}/snapshot_{snap_str}.hdf5')
print()
print("loaded")
print()
caes_obj = caesar.CAESAR(ds)
print()
print("object setup")
print()
#caes_obj.member_search(haloid = 'fof',blackholes=True,skipran=False,fof6d = True,nproc = pro)
caes_obj.member_search(haloid='fof', fof6d_file=fof_save_loc, fsps_bands='uvoir',
        ssp_model='FSPS', ssp_table_file=ssp_table_file,
        ext_law='composite', nproc=pro)
print()
print("member search")
print()
caes_obj.save(caesar_save_loc)

print()
print("saved")
print()



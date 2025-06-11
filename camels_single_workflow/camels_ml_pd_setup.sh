#!/bin/bash
#SBATCH --job-name=camels_setup
#SBATCH --output=output.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=d.zimmerman@ufl.edu
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=60gb
#SBATCH --account=narayanan
#SBATCH --qos=narayanan
#SBATCH --time=96:00:00
#SBATCH --array=14,18,24,32,44,62,90



date;hostname;pwd;
cd /home/d.zimmerman
module purge
source .bashrc

conda activate /blue/narayanan/d.zimmerman/code_environments/master_el9_environment

module load git
#module load gcc/12.2.0
module load intel/2025.1.0
module load openmpi/5.0.7
module load hdf5/1.14.6


sim=$(awk '{if(NR==1) print $0}' /orange/narayanan/d.zimmerman/camels_scripts/camels_config.txt)
run_1p_param=$(awk '{if(NR==2) print $0}' /orange/narayanan/d.zimmerman/camels_scripts/camels_config.txt)
run_1p_num=$(awk '{if(NR==3) print $0}' /orange/narayanan/d.zimmerman/camels_scripts/camels_config.txt)

caes_in_name="caesar_newsnaps"
caes_out_name="caesar"
subfind_in_name="groups"
subfind_out_name="groups"

if [ "$sim" = "0" ]
then
        sim_name="SIMBA"
	mkdir /orange/narayanan/d.zimmerman/camels_results/filtered/${sim_name}_1P_p${run_1p_param}_${run_1p_num}/snap${SLURM_ARRAY_TASK_ID}
elif [ "$sim" = "1" ]
then
        sim_name="IllustrisTNG"
	mkdir /orange/narayanan/d.zimmerman/camels_results/filtered/${sim_name}_1P_p${run_1p_param}_${run_1p_num}/snap${SLURM_ARRAY_TASK_ID}
else
        sim_name="failed - 1 or 0 not entered"
fi

echo "SIM TYPE:$sim_name"

sim_id="${sim_name}_1P_p${run_1p_param}_${run_1p_num}"

echo $SLURM_ARRAY_TASK_ID

if [ $(($SLURM_ARRAY_TASK_ID)) -lt 10 ]
then
	cp /orange/narayanan/d.zimmerman/camels_results/catalogs_loaded/${sim_id}/${caes_in_name}_00${SLURM_ARRAY_TASK_ID}.hdf5 /orange/narayanan/d.zimmerman/camels_results/catalogs_saved/${sim_id}/${caes_out_name}_${SLURM_ARRAY_TASK_ID}.hdf5
	cp /orange/narayanan/d.zimmerman/camels_results/catalogs_loaded/${sim_id}/${subfind_in_name}_00${SLURM_ARRAY_TASK_ID}.hdf5 /orange/narayanan/d.zimmerman/camels_results/catalogs_saved/${sim_id}/${subfind_out_name}_${SLURM_ARRAY_TASK_ID}.hdf5
elif [ $(($SLURM_ARRAY_TASK_ID)) -lt 100 ]
then
	cp /orange/narayanan/d.zimmerman/camels_results/catalogs_loaded/${sim_id}/${caes_in_name}_0${SLURM_ARRAY_TASK_ID}.hdf5 /orange/narayanan/d.zimmerman/camels_results/catalogs_saved/${sim_id}/${caes_out_name}_${SLURM_ARRAY_TASK_ID}.hdf5
	cp /orange/narayanan/d.zimmerman/camels_results/catalogs_loaded/${sim_id}/${subfind_in_name}_0${SLURM_ARRAY_TASK_ID}.hdf5 /orange/narayanan/d.zimmerman/camels_results/catalogs_saved/${sim_id}/${subfind_out_name}_${SLURM_ARRAY_TASK_ID}.hdf5
else
        cp /orange/narayanan/d.zimmerman/camels_results/catalogs_loaded/${sim_id}/${caes_in_name}_${SLURM_ARRAY_TASK_ID}.hdf5 /orange/narayanan/d.zimmerman/camels_results/catalogs_saved/${sim_id}/${caes_out_name}_${SLURM_ARRAY_TASK_ID}.hdf5
	cp /orange/narayanan/d.zimmerman/camels_results/catalogs_loaded/${sim_id}/${subfind_in_name}_${SLURM_ARRAY_TASK_ID}.hdf5 /orange/narayanan/d.zimmerman/camels_results/catalogs_saved/${sim_id}/${subfind_out_name}_${SLURM_ARRAY_TASK_ID}.hdf5
fi


if [ "${sim_name}" = "SIMBA" ]
then
	python /orange/narayanan/d.zimmerman/camels_scripts/filter_simba_camels_setup.py $sim_name $run_1p_param $run_1p_num $SLURM_ARRAY_TASK_ID -1
	python /orange/narayanan/d.zimmerman/camels_scripts/galaxy_positions.py ${sim_id} ${SLURM_ARRAY_TASK_ID}
	python /orange/narayanan/d.zimmerman/camels_scripts/caesar_good_gal_script.py ${sim_id} ${SLURM_ARRAY_TASK_ID}
	python /orange/narayanan/d.zimmerman/camels_scripts/powderday_setup.py ${sim_id} ${SLURM_ARRAY_TASK_ID}
	cp /orange/narayanan/d.zimmerman/camels_scripts/parameters_master_camels_simba.py /orange/narayanan/d.zimmerman/camels_results/pd_scripts/${sim_id}/snap${SLURM_ARRAY_TASK_ID}/parameters_master_camels.py
else
	python /orange/narayanan/d.zimmerman/camels_scripts/filter_tng_camels_setup.py $sim_name $run_1p_param $run_1p_num $SLURM_ARRAY_TASK_ID -1
	python /orange/narayanan/d.zimmerman/camels_scripts/galaxy_positions.py ${sim_id} ${SLURM_ARRAY_TASK_ID}
	python /orange/narayanan/d.zimmerman/camels_scripts/subfind_good_gal_script.py ${sim_id} ${SLURM_ARRAY_TASK_ID}
	python /orange/narayanan/d.zimmerman/camels_scripts/powderday_setup.py ${sim_id} ${SLURM_ARRAY_TASK_ID}
	cp /orange/narayanan/d.zimmerman/camels_scripts/parameters_master_camels_illustris.py /orange/narayanan/d.zimmerman/camels_results/pd_scripts/${sim_id}/snap${SLURM_ARRAY_TASK_ID}/parameters_master_camels.py
fi




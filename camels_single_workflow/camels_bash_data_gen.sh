#!/bin/bash
#SBATCH --job-name=camels_ml_data_gen
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

#V0 bash simulation loader/requester for CAMELS data
#this assumes you have made a Globus guest collection on hipergator to run things on
#SIMBA or IllustrisTNG 0 or 1

#camels globus ID

date;hostname;pwd;
cd /home/d.zimmerman
module purge
source .bashrc

conda activate /blue/narayanan/d.zimmerman/code_environments/master_el8_environment

module load git
#module load gcc/12.2.0
module load intel/2020.0.166
module load openmpi/4.1.5
module load hdf5/1.14.1

#Notes of interest:
LOCAL_SIMS_DIR="/orange/narayanan/d.zimmerman/camels_results/sims_loaded/"
LOCAL_CAESAR_DIR="/orange/narayanan/d.zimmerman/camels_results/catalogs_loaded/"
sim_type=$(awk '{if(NR==1) print $0}' /orange/narayanan/d.zimmerman/camels_scripts/camels_config.txt)
run_1p_param=$(awk '{if(NR==2) print $0}' /orange/narayanan/d.zimmerman/camels_scripts/camels_config.txt)
run_1p_num=$(awk '{if(NR==3) print $0}' /orange/narayanan/d.zimmerman/camels_scripts/camels_config.txt)

echo "$sim_type"

if [ "$sim_type" = "0" ]
then
	sim_name="SIMBA"
elif [ "$sim_type" = "1" ]
then
    	sim_name="IllustrisTNG"
else
	sim_name="failed - 1 or 0 not entered"
fi

echo "SIM TYPE:$sim_name"

sim_id="${sim_name}_1P_p${run_1p_param}_${run_1p_num}"

sim_dir="${LOCAL_SIMS_DIR}${sim_id}/"

cd /home/d.zimmerman/ml_sed_fitter/
python dataset_from_sim_gen.py $sim_id 67.11 $SLURM_ARRAY_TASK_ID


old_caes_dir="${LOCAL_CAESAR_DIR}${sim_id}/"


cd $LOCAL_SIMS_DIR
rm -r ${sim_id}/*
rmdir ${sim_id}


cd $LOCAL_CAESAR_DIR
rm -r ${sim_id}/*
rmdir ${sim_id}

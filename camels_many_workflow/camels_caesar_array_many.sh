#!/bin/bash
#SBATCH --job-name=camels_caesar_gen
#SBATCH --output=array_caesar_%A-%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=d.zimmerman@ufl.edu
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=60gb
#SBATCH --account=paul.torrey
#SBATCH --qos=paul.torrey-b
#SBATCH --time=10:00:00
#SBATCH --array=14,18,24,32,44,62,90

date;hostname;pwd;
cd /home/d.zimmerman
module purge
source .bashrc

conda activate /blue/narayanan/d.zimmerman/code_environments/master_el9_environment

module load git
module load intel/2025.1.0
module load openmpi/5.0.7
module load hdf5/1.14.6

cd /orange/narayanan/d.zimmerman/camels_scripts/

sim_num=$1
sim=$(awk '{if(NR==4*n+1) print $0}' n=${sim_num} /orange/narayanan/d.zimmerman/camels_scripts/camels_config2.txt)
run_1p_param=$(awk '{if(NR==4*n+2) print $0}' n=${sim_num} /orange/narayanan/d.zimmerman/camels_scripts/camels_config2.txt)
run_1p_num=$(awk '{if(NR==4*n+3) print $0}' n=${sim_num} /orange/narayanan/d.zimmerman/camels_scripts/camels_config2.txt)


if [ "$sim" = "0" ]
then
        sim_name="SIMBA"
elif [ "$sim" = "1" ]
then
        sim_name="IllustrisTNG"
else
        sim_name="failed - 1 or 0 not entered"
fi

echo "SIM TYPE:$sim_name"

sim_id="${sim_name}_1P_p${run_1p_param}_${run_1p_num}"

echo "SIM ID:$sim_id"


python camels_caesar_gen.py ${sim_name} ${run_1p_param} ${run_1p_num} ${SLURM_ARRAY_TASK_ID}

date


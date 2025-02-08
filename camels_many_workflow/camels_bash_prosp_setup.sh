#!/bin/bash
cd /home/d.zimmerman/
source .bashrc
module purge
conda deactivate
conda activate /blue/narayanan/d.zimmerman/code_environments/master_el8_environment
module load git
module load intel/2020.0.166
module load openmpi/4.1.5
module load hdf5/1.14.1
ml
cd /orange/narayanan/d.zimmerman/camels_scripts/
snum_min=$1
snum_max=$2


for ((cur_snum=${snum_min}; cur_snum<=${snum_max}; cur_snum++))
do
	sim_num=${cur_snum}
	sim=$(awk '{if(NR==4*n+1) print $0}' n=${sim_num} /orange/narayanan/d.zimmerman/camels_scripts/camels_config2.txt)
	run_1p_param=$(awk '{if(NR==4*n+2) print $0}' n=${sim_num} /orange/narayanan/d.zimmerman/camels_scripts/camels_config2.txt)
	run_1p_num=$(awk '{if(NR==4*n+3) print $0}' n=${sim_num} /orange/narayanan/d.zimmerman/camels_scripts/camels_config2.txt)

	if [ "$sim" = "0" ]
	then
		sim_name="SIMBA"
	else [ "$sim" = "1" ]
		sim_name="IllustrisTNG"
	fi
	
	sim_id="${sim_name}_1P_p${run_1p_param}_${run_1p_num}"
	echo "SIM ID:$sim_id"
	# 14 18 24 32 44 62 90
	for snap in 14
	do
		echo ${snap}
		python "/orange/narayanan/d.zimmerman/camels_scripts/prospector_setup.py" ${sim_id} ${snap}
		# ${sim_id} ${snap}"
	done
	echo -e "\n"
	
done	

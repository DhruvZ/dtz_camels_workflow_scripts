#!/bin/bash

sim=$(awk '{if(NR==1) print $0}' /orange/narayanan/d.zimmerman/camels_scripts/camels_config.txt)
run_1p_param=$(awk '{if(NR==2) print $0}' /orange/narayanan/d.zimmerman/camels_scripts/camels_config.txt)
run_1p_num=$(awk '{if(NR==3) print $0}' /orange/narayanan/d.zimmerman/camels_scripts/camels_config.txt)

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

# snap list: 0,1,2,4,10,19,33
# 14 18 24 32 44 62 90


for snum in 62 90 
do
	echo "on snap ${snum}"
	path=/orange/narayanan/d.zimmerman/camels_results/pd_scripts/${sim_id}/snap${snum}/master.snap${snum}.job
	echo $path
	cd /orange/narayanan/d.zimmerman/camels_results/pd_scripts/${sim_id}/snap${snum}/
	sbatch $path
	echo -e "\n"

done
#sbatch 


#!/bin/bash 

#V1 bash simulation loader/requester for CAMELS data
#this assumes you have made a Globus guest collection on hipergator to run things on
#SIMBA or IllustrisTNG 0 or 1

#camels globus ID

#Notes of interest:
CAMELS_GLOBUS="58bdcd24-6590-11ec-9b60-f9dfb1abb183"
HPG_COLLECT_GLOBUS="105bd89d-82b6-4957-be97-9bd47f628101"
LOCAL_SIMS_DIR="/sims_loaded/"
LOCAL_CAESAR_DIR="/catalogs_loaded/"
CAMELS_SIMS_DIR="/Sims/"
CAMELS_CAESAR_DIR="/Caesar/"
CAMELS_SUBFIND_DIR="/FOF_Subfind/"

snum_min=$1
snum_max=$2
caesar_load=$3

module load globus

for ((cur_snum=${snum_min}; cur_snum<=${snum_max}; cur_snum++))
do
	sim_type=$(awk '{if(NR==4*n+1) print $0}' n=${cur_snum} /orange/narayanan/d.zimmerman/camels_scripts/camels_config2.txt)
	run_1p_param=$(awk '{if(NR==4*n+2) print $0}' n=${cur_snum} /orange/narayanan/d.zimmerman/camels_scripts/camels_config2.txt)
	run_1p_num=$(awk '{if(NR==4*n+3) print $0}' n=${cur_snum} /orange/narayanan/d.zimmerman/camels_scripts/camels_config2.txt)
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
	
	in_sim_dir="${CAMELS_SIMS_DIR}${sim_name}/1P/1P_p${run_1p_param}_${run_1p_num}/"
	in_caes_dir="${CAMELS_CAESAR_DIR}${sim_name}/L25n256/1P/1P_p${run_1p_param}_${run_1p_num}/"
	in_subfind_dir="${CAMELS_SUBFIND_DIR}${sim_name}/1P/1P_p${run_1p_param}_${run_1p_num}/"
	out_sim_dir="/sims_loaded/${sim_name}_1P_p${run_1p_param}_${run_1p_num}/"
	out_cat_dir="/catalogs_loaded/${sim_name}_1P_p${run_1p_param}_${run_1p_num}/"
	
	
	globus transfer --recursive --dry-run ${CAMELS_GLOBUS}:${in_sim_dir} ${HPG_COLLECT_GLOBUS}:${out_sim_dir}
	echo -e "\n"
	globus transfer --recursive --dry-run ${CAMELS_GLOBUS}:${in_caes_dir} ${HPG_COLLECT_GLOBUS}:${out_cat_dir}
	echo -e "\n"
	globus transfer --recursive --dry-run ${CAMELS_GLOBUS}:${in_subfind_dir} ${HPG_COLLECT_GLOBUS}:${out_cat_dir}
	echo -e "\n"


	echo -e "Please confirm the transfer is correct (y/n):"
	
	read confirm
	
	if [ "$confirm" = "y" ]
	then
		echo -e "confirmed"
		globus transfer --recursive ${CAMELS_GLOBUS}:${in_sim_dir} ${HPG_COLLECT_GLOBUS}:${out_sim_dir}
		echo -e "\n"
		if [ "$caesar_load" = "1" ]
		then
			globus transfer --recursive ${CAMELS_GLOBUS}:${in_caes_dir} ${HPG_COLLECT_GLOBUS}:${out_cat_dir}
			echo -e "loaded caesar files"
			echo -e "\n"
		else
			echo -e "no caesar load"
			echo -e "\n"
		fi
				
		globus transfer --recursive ${CAMELS_GLOBUS}:${in_subfind_dir} ${HPG_COLLECT_GLOBUS}:${out_cat_dir}
		echo -e "\n"
		
		mkdir /orange/narayanan/d.zimmerman/camels_results/filtered/${sim_name}_1P_p${run_1p_param}_${run_1p_num}/
		mkdir /orange/narayanan/d.zimmerman/camels_results/pd_scripts/${sim_name}_1P_p${run_1p_param}_${run_1p_num}/
		mkdir /orange/narayanan/d.zimmerman/camels_results/pd_runs/${sim_name}_1P_p${run_1p_param}_${run_1p_num}/
		mkdir /orange/narayanan/d.zimmerman/camels_results/catalogs_saved/${sim_name}_1P_p${run_1p_param}_${run_1p_num}/
	else
		echo -e "transfer cancelled"
	fi
done	

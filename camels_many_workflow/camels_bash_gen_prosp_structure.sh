snum_min=$1
snum_max=$2


base_path="/orange/narayanan/d.zimmerman/camels_results/prosp_runs/"

for ((cur_snum=${snum_min}; cur_snum<=${snum_max}; cur_snum++))
do
        sim_type=$(awk '{if(NR==4*n+1) print $0}' n=${cur_snum} /orange/narayanan/d.zimmerman/camels_scripts/camels_config2.txt)
        run_1p_param=$(awk '{if(NR==4*n+2) print $0}' n=${cur_snum} /orange/narayanan/d.zimmerman/camels_scripts/camels_config2.txt)
        run_1p_num=$(awk '{if(NR==4*n+3) print $0}' n=${cur_snum} /orange/narayanan/d.zimmerman/camels_scripts/camels_config2.txt)

        if [ "$sim_type" = "0" ]
        then
                sim_name="SIMBA"
        elif [ "$sim_type" = "1" ]
        then
                sim_name="IllustrisTNG"
        else
                sim_name="failed - 1 or 0 not entered"
        fi

	sim_id="${sim_name}_1P_p${run_1p_param}_${run_1p_num}"
	echo -e "\n"
	echo "SIM ID:$sim_id"

	mkdir "${base_path}${sim_id}"

	for snap in 14 18 24 32 44 62 90
	do
		echo -e "${base_path}${sim_id}/snap${snap}"
		echo -e "${base_path}${sim_id}/snap${snap}/all_photometry/"
		echo -e "${base_path}${sim_id}/snap${snap}/limited_photometry/"
		echo -e "\n"
		mkdir "${base_path}${sim_id}/snap${snap}"
                mkdir "${base_path}${sim_id}/snap${snap}/all_photometry/"
                mkdir "${base_path}${sim_id}/snap${snap}/limited_photometry/"

	done
done

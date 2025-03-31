snum_min=$1
snum_max=$2


for ((cur_snum=${snum_min}; cur_snum<=${snum_max}; cur_snum++))
do
	echo "on number:${cur_snum}"
	sbatch camels_many_data_gen.sh ${cur_snum}
done	

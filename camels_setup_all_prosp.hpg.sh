#!/bin/bash 

#PROSPECTOR cluster setup convenience script for SLURM queue manager
#on HiPerGator at the University of FLorida.  This sets up the model
#files for a cosmological simulation where we want to model many
#galaxies at once.

#> echo $BASH_VERSION

base_path=${1}
sim_id=${2}
snap=${3}
max_gal_id=${4}

echo "processing model file for sim,snapshot: $sim_id,$snap"


#clear the pyc files
rm -f *.pyc

job_dir="${base_path}/prosp_runs/${sim_id}/snap${snap}/"

echo "writing slurm submission script file"
qsubfile="${job_dir}prosp_master.${sim_id}_snap${snap}.job"
rm -f $qsubfile
echo $qsubfile

echo "#! /bin/bash" >>$qsubfile
echo "#SBATCH --job-name=${sim_id}.snap${snap}_prosp" >>$qsubfile
echo "#SBATCH --output=${sim_id}.snap${snap}_%A-%a.o" >>$qsubfile
echo "#SBATCH --error=${sim_id}.snap${snap}_%A-%a.e" >>$qsubfile
echo "#SBATCH --mail-type=ALL" >>$qsubfile
echo "#SBATCH --mail-user=d.zimmerman@ufl.edu" >>$qsubfile
echo "#SBATCH --time=48:00:00" >>$qsubfile
echo "#SBATCH --mem=10gb">>$qsubfile
echo "#SBATCH --account=narayanan">>$qsubfile
echo "#SBATCH --qos=narayanan-b">>$qsubfile
echo "#SBATCH --array=0-${max_gal_id}">>$qsubfile
echo -e "\n">>$qsubfile
echo -e "\n" >>$qsubfile

echo "cd /home/d.zimmerman">>$qsubfile
echo "module purge">>$qsubfile
echo "source .bashrc">>$qsubfile
echo "conda activate /blue/narayanan/d.zimmerman/code_environments/master_el8_environment">>$qsubfile
echo -e "\n">>$qsubfile
echo "module load git">>$qsubfile
echo "module load intel/2020.0.166">>$qsubfile
echo "module load openmpi/4.1.6">>$qsubfile
echo "module load hdf5/1.14.1">>$qsubfile
echo -e "\n">>$qsubfile


echo "sim_id=\"${sim_id}\"">>$qsubfile
echo "snap=\"${snap}\"">>$qsubfile
echo "gal_num=\${SLURM_ARRAY_TASK_ID}">>$qsubfile
echo "full_phot=\${1}">>$qsubfile
echo "z_idx=\${2}">>$qsubfile
echo "filt_file=\${3}">>$qsubfile
echo "snr=\${4}">>$qsubfile
echo "extra_label=\${5}">>$qsubfile
echo "root_override=\"${base_path}\"">>$qsubfile
echo "z_known=\${6}">>$qsubfile

echo "cd \"/orange/narayanan/d.zimmerman/camels_scripts/\"">>$qsubfile

echo -e "if [ \"\$z_known\" = \"1\" ]">>$qsubfile
echo -e "then">>$qsubfile
echo -e "\tpython run_prosp_ml_sed.py \${sim_id} \${snap} \${gal_num} \${full_phot} \${z_idx} \${filt_file} \${snr} \${extra_label} \${root_override}">>$qsubfile
echo -e "\tpython process_prosp_output_camels_dtz.py \${sim_id} \${snap} \${gal_num} \${full_phot} \${z_idx} \${extra_label} \${root_override}">>$qsubfile
echo -e "else">>$qsubfile
echo -e "\tpython run_prosp_ml_sed_zvar.py \${sim_id} \${snap} \${gal_num} \${full_phot} \${z_idx} \${filt_file} \${snr} \${extra_label} \${root_override}">>$qsubfile
echo -e "\tpython process_prosp_output_camels_zvar_dtz.py \${sim_id} \${snap} \${gal_num} \${full_phot} \${z_idx} \${extra_label} \${root_override}">>$qsubfile
echo -e "fi">>$qsubfile

echo "date">>$qsubfile
#done


#!/bin/bash
#SBATCH --job-name=slurm_test
#SBATCH --output=stest.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=d.zimmerman@ufl.edu
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=60gb
#SBATCH --account=narayanan
#SBATCH --qos=narayanan
#SBATCH --time=96:00:00


pwd
echo $1
echo $2

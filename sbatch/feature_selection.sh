#!/bin/sh
#SBATCH --job-name=feature_selection
#SBATCH --output=/gpfs/data/chiu-lab/bw_folder/dissertation/5hmc_tx/interation_100_allrace_newMM/logs/res_%A.out 
#SBATCH --error=/gpfs/data/chiu-lab/bw_folder/dissertation/5hmc_tx/interation_100_allrace_newMM/logs/res_%A.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=23gb
#SBATCH --partition=tier1q
i=${ARGS1}
module load gcc/12.1.0
module load miniconda3/23.1.0
source activate /home/beiw/miniconda3/envs/5hmc_os_env
echo "Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"

Rscript /gpfs/data/chiu-lab/bw_folder/dissertation/5hmc_tx/interation_100_allrace_newMM/r_scripts/feature_selection.R ${i}
#!/bin/bash -x
#SBATCH -A naiss2024-5-311
#SBATCH --output=/proj/elofssonlab/users/x_safro/logs/%A_%a.out
#SBATCH --error=/proj/elofssonlab/users/x_safro/logs/%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --export=ALL,CUDA_VISIBLE_DEVICES

# Load AlphaFold module
module load Mambaforge/23.3.1-1-hpc1

conda activate af_data_env

python rename_chai1_files.py

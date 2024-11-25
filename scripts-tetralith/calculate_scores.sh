#!/bin/bash -x
#SBATCH -A naiss2024-5-311
#SBATCH --output=/proj/elofssonlab/users/x_safro/logs/%A_%a.out
#SBATCH --error=/proj/elofssonlab/users/x_safro/logs/%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=32
#SBATCH -t 08:00:00
#SBATCH --export=ALL,CUDA_VISIBLE_DEVICES

# Load AlphaFold module
module load Miniforge/24.7.1-2-hpc1 

mamba activate snakemake_env

BASEDIR="/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/"
CONFIGFILE="$BASEDIR/benchmarking-workflow/config/config_prebuilt.yaml"


snakemake -s calculate_scores.smk --use-conda --configfile config/config_prebuilt.yaml --directory $BASEDIR
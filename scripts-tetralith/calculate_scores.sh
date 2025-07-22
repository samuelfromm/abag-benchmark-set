#!/bin/bash -x
#SBATCH -A naiss2024-5-311
#SBATCH --output=/proj/elofssonlab/users/x_safro/logs/%A_%a.out
#SBATCH --error=/proj/elofssonlab/users/x_safro/logs/%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=32
#SBATCH -t 48:00:00
#SBATCH --export=ALL,CUDA_VISIBLE_DEVICES

# Load AlphaFold module
module load Miniforge/24.7.1-2-hpc1 

mamba activate snakemake_env



BASEDIR="/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/benchmarking-workflow"
RUNNAME="default"
CONFIGFILE="$BASEDIR/config/config_prebuilt.yaml"
OUTDIR="/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/scores/$RUNNAME/output"
SAMPLES_CSV="/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/scores/$RUNNAME/input/samples_$RUNNAME.csv"
SMKFILE="$BASEDIR/calculate_scores.smk"
MERGESMKFILE="$BASEDIR/merge.smk"

#snakemake -s $SMKFILE --use-conda --configfile $CONFIGFILE --directory $BASEDIR --cores $SLURM_CPUS_PER_TASK --config output_dir=$OUTDIR samples_csv=${SAMPLES_CSV} --keep-going

snakemake -s $MERGESMKFILE --use-conda --configfile $CONFIGFILE --directory $BASEDIR --cores $SLURM_CPUS_PER_TASK --config output_dir=$OUTDIR samples_csv=${SAMPLES_CSV} --keep-going
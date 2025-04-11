#!/bin/bash -x
#SBATCH -A naiss2024-5-311
#SBATCH --output=/proj/elofssonlab/users/x_safro/logs/%A_%a.out
#SBATCH --error=/proj/elofssonlab/users/x_safro/logs/%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=8
#SBATCH -t 12:00:00
#SBATCH --export=ALL,CUDA_VISIBLE_DEVICES

module load Miniforge/24.7.1-2-hpc1

mamba activate snakemake_env

# File containing IDs to process, one per line
IDFILE="/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/IDs_scores_missing.csv"

# Ensure the ID file exists
if [[ ! -f "$IDFILE" ]]; then
  echo "Error: ID file $IDFILE not found!"
  exit 1
fi

# Fetch the PDBID corresponding to this SLURM array task
PDBID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$IDFILE")


# Check if PDBID was successfully fetched
if [[ -z "$PDBID" ]]; then
  echo "Error: No PDBID found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

RUNNAME="default"

BASEDIR="/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/benchmarking-workflow"
CONFIGFILE="$BASEDIR/config/config_actifptm_prebuilt.yaml"
OUTDIR="/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/scores_actifptm/$RUNNAME/output/$PDBID"
SAMPLES_CSV="/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/scores_actifptm/$RUNNAME/input/${PDBID}_samples_$RUNNAME.csv"
SMKFILE="$BASEDIR/calculate_scores_actifptm.smk"

# Ensure necessary directories exist
mkdir -p "$OUTDIR"

# Ensure the samples CSV file exists
if [[ ! -f "$SAMPLES_CSV" ]]; then
  echo "Error: Samples CSV file $SAMPLES_CSV not found for PDBID=$PDBID!"
  exit 1
fi

# Run Snakemake
snakemake -s "$SMKFILE" \
  --use-conda \
  --configfile "$CONFIGFILE" \
  --directory "$BASEDIR" \
  --cores "$SLURM_CPUS_PER_TASK" \
  --config output_dir="$OUTDIR" samples_csv="$SAMPLES_CSV" \
  --keep-going
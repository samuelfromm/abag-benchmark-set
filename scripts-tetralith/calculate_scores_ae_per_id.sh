#!/bin/bash -x
#SBATCH -A Berzelius-2025-214
#SBATCH --output=/proj/berzelius-2021-29/users/x_safro/git/paper/logs/%A_%a.out
#SBATCH --error=/proj/berzelius-2021-29/users/x_safro/git/paper/logs/%A_%a.err
#SBATCH --array=1-9
#SBATCH --cpus-per-task=8
#SBATCH --partition=berzelius-cpu
#SBATCH -t 12:00:00
#SBATCH --export=ALL,CUDA_VISIBLE_DEVICES


module load Mambaforge/23.3.1-1-hpc1-bdist

mamba activate snakemake_env


# File containing IDs to process, one per line
IDFILE="/proj/berzelius-2021-29/users/x_safro/git/paper/abag-benchmark-set/IDs_scores_missing.csv"

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

RUNNAME="alphafold3"

BASEDIR="/proj/berzelius-2021-29/users/x_safro/git/paper/abag-benchmark-set/benchmarking-workflow"
CONFIGFILE="$BASEDIR/config/config_ae_prebuilt.yaml"
OUTDIR="/proj/berzelius-2021-29/users/x_safro/git/paper/abag-benchmark-set/data/scores_ae/$RUNNAME/output/$PDBID"
SAMPLES_CSV="/proj/berzelius-2021-29/users/x_safro/git/paper/abag-benchmark-set/data/scores_ae/$RUNNAME/input/${PDBID}_samples_$RUNNAME.csv"
SMKFILE="$BASEDIR/calculate_scores_aligned_error.smk"

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
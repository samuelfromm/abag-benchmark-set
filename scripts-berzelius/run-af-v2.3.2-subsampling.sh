#!/bin/bash -x
#SBATCH -A berzelius-2024-220
#SBATCH --output=/proj/berzelius-2021-29/users/x_safro/git/abag-benchmark-set/logs/%A_%a.out
#SBATCH --error=/proj/berzelius-2021-29/users/x_safro/git/abag-benchmark-set/logs/%A_%a.err
#SBATCH --array=1
#SBATCH --gpus=1
#SBATCH -t 60:00:00
#SBATCH --export=ALL,CUDA_VISIBLE_DEVICES

# Load necessary modules
module load AlphaFold/2.3.2-hpc1

# Set offset based on input parameter or default to zero
offset=${1:-0}
LN=$(( SLURM_ARRAY_TASK_ID + offset ))

RUNNAME="subsampling-16-32"

# File paths
IDFILE="/proj/berzelius-2021-29/users/x_safro/git/abag-benchmark-set/run-info/IDs_${RUNNAME}_missing.csv"
PDBID=$(sed -n "${LN}p" $IDFILE)


FASTAPATH="/proj/berzelius-2021-29/users/x_safro/git/abag-benchmark-set/data/db/complexFastas/${PDBID}.fasta"
OUTPUTDIR="/proj/berzelius-2021-29/users/x_safro/git/abag-benchmark-set/data/models/${RUNNAME}"
MSADIR="/proj/berzelius-2021-29/users/x_safro/git/abag-benchmark-set/data/MSAs"

# Set directories for AlphaFold
AFHOME="/proj/berzelius-2021-29/users/x_safro/programs/alphafold-clone/"
DATABASE_DIR="/proj/berzelius-2021-29/alphafold_data_v2.3"
PARAMS_DIR=$DATABASE_DIR

# AlphaFold configuration flags
flags=(
  "--uniprot_database_path=${DATABASE_DIR}/uniprot/uniprot.fasta"
  "--uniref90_database_path=${DATABASE_DIR}/uniref90/uniref90.fasta"
  "--mgnify_database_path=${DATABASE_DIR}/mgnify/mgy_clusters_2022_05.fa"
  "--bfd_database_path=${DATABASE_DIR}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
  "--pdb_seqres_database_path=${DATABASE_DIR}/pdb_seqres/pdb_seqres.txt"
  "--uniref30_database_path=${DATABASE_DIR}/uniref30/UniRef30_2021_03"
  "--template_mmcif_dir=${DATABASE_DIR}/pdb_mmcif/mmcif_files"
  "--obsolete_pdbs_path=${DATABASE_DIR}/pdb_mmcif/obsolete.dat"
  "--db_preset=full_dbs"
  "--use_precomputed_msas=True"
  "--fasta_paths=${FASTAPATH}"
  "--output_dir=${OUTPUTDIR}"
  "--msa_dir=${MSADIR}"
  "--data_dir=${PARAMS_DIR}"
  "--max_template_date=2021-09-30"
  "--use_gpu_relax=False"
  "--run_features_only=False"
  "--rank_models=False"
)

run_specific_flags=(
  "--model_preset=multimer" 
  "--models_to_relax=none"
  "--num_multimer_predictions_per_model=40"
  "--starting_prediction_per_model=0"
  "--suffix=${RUNNAME}"
  "--reuse_features=True"
  "--num_msa=16"
  "--num_extra_msa=32"
  "--resample_msa_in_recycling=True"
)


# Start timer
SECONDS=0

# Run AlphaFold
python3 "${AFHOME}/run_alphafold.py" "${flags[@]}" "${run_specific_flags[@]}"

# Log elapsed time
echo "ELAPSED TIME: $(($SECONDS / 3600))h $((($SECONDS / 60) % 60))m $(($SECONDS % 60))s"
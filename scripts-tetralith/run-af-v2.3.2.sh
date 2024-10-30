#!/bin/bash -x
#SBATCH -A naiss2024-5-311
#SBATCH --output=/proj/elofssonlab/users/x_safro/logs/%A_%a.out
#SBATCH --error=/proj/elofssonlab/users/x_safro/logs/%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=12
#SBATCH -t 06:00:00
#SBATCH --export=ALL,CUDA_VISIBLE_DEVICES


module load AlphaFold/2.3.2-hpc1


if [ -z $1 ]
then
        offset=0
else
        offset=$1
fi


LN=$(( SLURM_ARRAY_TASK_ID + offset ))


IDFILE=/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/db/IDs.csv
PDBID=$(sed -n ${LN}p $IDFILE) 

#PDBID="2ORB"



FASTAPATH=/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/db/complexFastas/${PDBID}complex.fasta
OUTPUTDIR=/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/MSAs


AFHOME=/proj/elofssonlab/users/x_safro/programs/alphafold-clone/

DATABASE_DIR=/proj/elofssonlab/alphafold_data_v2.3
PARAMS_DIR=/proj/elofssonlab/alphafold_data_v2.3


flags="
--uniprot_database_path=$DATABASE_DIR/uniprot/uniprot.fasta
--uniref90_database_path=$DATABASE_DIR/uniref90/uniref90.fasta
--mgnify_database_path=$DATABASE_DIR/mgnify/mgy_clusters_2022_05.fa 
--bfd_database_path=$DATABASE_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
--pdb_seqres_database_path=$DATABASE_DIR/pdb_seqres/pdb_seqres.txt
--uniref30_database_path=$DATABASE_DIR/uniref30/UniRef30_2021_03
--template_mmcif_dir=$DATABASE_DIR/pdb_mmcif/mmcif_files
--obsolete_pdbs_path=$DATABASE_DIR/pdb_mmcif/obsolete.dat
--db_preset=full_dbs
--use_precomputed_msas=True
--fasta_paths=$FASTAPATH
--output_dir=$OUTPUTDIR
--data_dir=$PARAMS_DIR
--max_template_date=2021-09-30
--use_gpu_relax=False
--model_preset=multimer
"

extra_flags="
--run_features_only=True
"


SECONDS=0

python3 $AFHOME/run_alphafold.py $flags $extra_flags

echo "ELAPSED TIME: $(($SECONDS/3600))h $((($SECONDS/60)%60))m $(($SECONDS%60))s"


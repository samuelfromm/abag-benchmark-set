

# Scoring

module load Mambaforge/23.3.1-1-hpc1
conda activate snakemake_env

snakemake -s calculate_scores.smk --use-conda --configfile config/config.yaml
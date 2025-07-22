#!/bin/bash


if [ -z "$1" ]; then
    echo "Usage: $0 <PRESET>"
    exit 1
fi

PRESET=$1   

INPUT_CSV_PATH="/home/sfromm/git/abag-benchmark-set/data/scores_ae/$PRESET/scores_$PRESET.csv"        # Path to the input CSV file
OUTPUT_CSV_PATH="/home/sfromm/git/abag-benchmark-set/results/sampling_scores_ae/sampling_scores_$PRESET.csv"  # Path to save the output results
MAX_SAMPLE_SIZE=200                   # Maximum sample size
STEP_SIZE=1                          # Step size for incrementing sample sizes
ITERATIONS=50                        # Number of iterations for averaging
COLUMNS_TO_ANALYZE="abag_dockq ranking_confidence min_pmidockq TM_normalized_reference aetm aeitm ae_ranking_confidence max_per_chain_aeitm min_per_chain_aeitm mean_per_chain_aeitm paetm paeitm pae_ranking_confidence max_per_chain_paeitm min_per_chain_paeitm mean_per_chain_paeitm min_per_chain_iptm mean_per_chain_iptm"  # Columns to analyze (space-separated)
REFERENCE_COLUMN="abag_dockq"         # Reference column to use for values


# Run the Python script with the defined parameters
python analysis/sample.py --input_csv_paths "$INPUT_CSV_PATH" \
                      --output_csv_path "$OUTPUT_CSV_PATH" \
                      --max_sample_size "$MAX_SAMPLE_SIZE" \
                      --step_size "$STEP_SIZE" \
                      --iterations "$ITERATIONS" \
                      --columns_to_analyze $COLUMNS_TO_ANALYZE \
                      --reference_column "$REFERENCE_COLUMN" \
                      --preset "$PRESET" \
                      --replace
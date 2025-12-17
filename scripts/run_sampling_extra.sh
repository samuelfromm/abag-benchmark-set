#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <PRESET>"
    exit 1
fi

PRESET=$1   

INPUT_CSV_PATH="/home/sfromm/git/abag-benchmark-set/data/scores_extra/$PRESET/scores_$PRESET.csv"        # Path to the input CSV file
OUTPUT_CSV_PATH="/home/sfromm/git/abag-benchmark-set/results/sampling_scores_extra/sampling_scores_$PRESET.csv"  # Path to save the output results
MAX_SAMPLE_SIZE=200                   # Maximum sample size
STEP_SIZE=1                          # Step size for incrementing sample sizes
ITERATIONS=50                        # Number of iterations for averaging
COLUMNS_TO_ANALYZE="abag_dockq,ptm,iptm,ranking_confidence,max_per_chain_iptm,min_per_chain_iptm,mean_per_chain_iptm,abag_max_ipSAE,abag_max_ipSAE_d0chn,abag_max_ipSAE_d0dom,abag_max_ipTM_d0chn,abag_max_pDockQ,abag_max_pDockQ2,abag_max_LIS,abag_asym_1_Chn1,abag_asym_1_ipSAE,abag_asym_1_ipSAE_d0chn,abag_asym_1_ipSAE_d0dom,abag_asym_1_ipTM_d0chn,abag_asym_1_pDockQ,abag_asym_1_pDockQ2,abag_asym_1_LIS,abag_asym_2_Chn1,abag_asym_2_ipSAE,abag_asym_2_ipSAE_d0chn,abag_asym_2_ipSAE_d0dom,abag_asym_2_ipTM_d0chn,abag_asym_2_pDockQ,abag_asym_2_pDockQ2,abag_asym_2_LIS"  # Columns to analyze (space-separated)
COLUMNS_TO_ANALYZE=${COLUMNS_TO_ANALYZE//,/ }  # Replace commas with spaces
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
#!/bin/bash

# Check if presets are provided
if [ "$#" -lt 1 ]; then
    echo "Error: No presets provided. Usage: $0 <preset1> [preset2 ...]"
    exit 1
fi

# Read presets from command-line arguments
PRESETS=("$@")

# Define constants
BASE_INPUT_PATH="/home/sfromm/git/abag-benchmark-set/data/scores"
BASE_OUTPUT_PATH="/home/sfromm/git/abag-benchmark-set/results/sampling_scores"
MAX_SAMPLE_SIZE=200
STEP_SIZE=1
ITERATIONS=50
COLUMNS_TO_ANALYZE="abag_dockq ranking_confidence min_pmidockq TM_normalized_reference"
REFERENCE_COLUMN="abag_dockq"

# Initialize a variable to store all input file paths
INPUT_CSV_PATHS=()
PRESET_NAMES=""

# Loop through presets to construct input paths and concatenate preset names
for PRESET in "${PRESETS[@]}"; do
    echo "Processing preset: $PRESET"
    INPUT_FILE="${BASE_INPUT_PATH}/${PRESET}/scores_${PRESET}.csv"

    if [ -f "$INPUT_FILE" ]; then
        INPUT_CSV_PATHS+=("$INPUT_FILE")
        PRESET_NAMES+="${PRESET}_"
    else
        echo "Warning: File not found for preset '$PRESET': $INPUT_FILE"
    fi
done

# Remove trailing underscore from PRESET_NAMES
PRESET_NAMES=${PRESET_NAMES%_}

# Set output file path
OUTPUT_CSV_PATH="${BASE_OUTPUT_PATH}/sampling_scores_${PRESET_NAMES}.csv"

# Check if at least one input file exists
if [ ${#INPUT_CSV_PATHS[@]} -eq 0 ]; then
    echo "Error: No valid input files found for presets."
    exit 1
fi

# Join all input file paths with spaces for the Python script
INPUT_CSV_PATHS_JOINED=$(printf " %s" "${INPUT_CSV_PATHS[@]}")

# Run the Python script with the constructed arguments
python analysis/sample.py --input_csv_paths ${INPUT_CSV_PATHS_JOINED} \
                          --output_csv_path "$OUTPUT_CSV_PATH" \
                          --max_sample_size "$MAX_SAMPLE_SIZE" \
                          --step_size "$STEP_SIZE" \
                          --iterations "$ITERATIONS" \
                          --columns_to_analyze $COLUMNS_TO_ANALYZE \
                          --reference_column "$REFERENCE_COLUMN" \
                          --preset "$PRESET_NAMES"

echo "Processing completed. Results saved to $OUTPUT_CSV_PATH"
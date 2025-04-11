#!/bin/bash

# Check if RUNNAME argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <RUNNAME>"
    exit 1
fi
# Set variables
RUNNAME="$1"
BASEDIR="/proj/elofssonlab/users/x_safro/git"
OUTDIR="$BASEDIR/abag-benchmark-set/data/scores/$RUNNAME/input"
NUMPREDICTIONS=40
NUMSAMPLES=5
DB_FILE="$BASEDIR/abag-benchmark-set/data/db/lightDb.txt"
MODELDIR="$BASEDIR/abag-benchmark-set/data/models"
FILTERED_PDB_DIR="$BASEDIR/abag-benchmark-set/data/db/structures_filtered"

# Templates for paths
PDB_TEMPLATE="$MODELDIR/$RUNNAME/{ID}complex/model_seed_{PRED_NUM}_sample_{MODEL_NUM}_${RUNNAME}.pdb" 
FEATURES_TEMPLATE="$MODELDIR/$RUNNAME/{ID}complex/full_confidences_seed_{PRED_NUM}_sample_{MODEL_NUM}_${RUNNAME}.json" #"$MODELDIR/$RUNNAME/{ID}complex/{ID}complex_data.json"
DATA_TEMPLATE="$MODELDIR/$RUNNAME/{ID}complex/full_confidences_seed_{PRED_NUM}_sample_{MODEL_NUM}_${RUNNAME}.json"
SAMPLE_ID_TEMPLATE="{ID}_seed_{PRED_NUM}_sample_{MODEL_NUM}_${RUNNAME}"

# Function to substitute variables in templates
substitute_template() {
    local template="$1"
    local id="$2"
    local model_num="$3"
    local pred_num="$4"

    template="${template//\{ID\}/$id}"
    template="${template//\{MODEL_NUM\}/$model_num}"
    template="${template//\{PRED_NUM\}/$pred_num}"
    echo "$template"
}

# Read the database file, skipping the header
tail -n +2 "$DB_FILE" | while IFS=, read -r pdb_id chain_A chain_H chain_L; do


    # Set output CSV file for the current pdb_id
    output_file="$OUTDIR/${pdb_id}_samples_$RUNNAME.csv"

    # Initialize the output file
    > "$output_file"
    echo "sample_id,pdbid,Achain,Hchain,Lchain,reference_pdb,query_pdb,query_af_features,query_af_data,model,prediction,preset" >> "$output_file"

    # Process predictions for the current pdb_id
    for ((MODEL_NUM = 1; MODEL_NUM <= NUMSAMPLES; MODEL_NUM++)); do
        for ((PRED_NUM = 0; PRED_NUM < NUMPREDICTIONS; PRED_NUM++)); do

            # Generate dynamic values
            SAMPLE_ID=$(substitute_template "$SAMPLE_ID_TEMPLATE" "$pdb_id" "$MODEL_NUM" "$PRED_NUM")
            QUERY_PDB=$(substitute_template "$PDB_TEMPLATE" "$pdb_id" "$MODEL_NUM" "$PRED_NUM")
            QUERY_DATA=$(substitute_template "$DATA_TEMPLATE" "$pdb_id" "$MODEL_NUM" "$PRED_NUM")
            QUERY_FEATURES=$(substitute_template "$FEATURES_TEMPLATE" "$pdb_id" "$MODEL_NUM" "$PRED_NUM")
            REFERENCE_PDB="$FILTERED_PDB_DIR/${pdb_id}_filtered.pdb"
            MODEL="model_${MODEL_NUM}_${WEIGHT_VERSION}"

            # Validate file existence
            all_files_exist=true
            for file in "$QUERY_PDB" "$REFERENCE_PDB" "$QUERY_DATA" "$QUERY_FEATURES"; do
                if [ ! -f "$file" ]; then
                    echo "Warning: File $file not found for $pdb_id"
                    all_files_exist=false
                    break
                fi
            done

            # Skip processing if any file is missing
            if [ "$all_files_exist" = false ]; then
                continue
            fi

            # Append the processed data to the output CSV file
            echo "${SAMPLE_ID},${pdb_id},${chain_A},${chain_H},${chain_L},${REFERENCE_PDB},${QUERY_PDB},${QUERY_FEATURES},${QUERY_DATA},${MODEL},${PRED_NUM},${RUNNAME}" >> "$output_file"
        done
    done
done
#!/bin/bash

# Check if RUNNAME argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <RUNNAME>"
    exit 1
fi

# Set variables
RUNNAME="$1"
BASEDIR="/proj/elofssonlab/users/x_safro/git"
OUTDIR="$BASEDIR/abag-benchmark-set/data/scores_paired/$RUNNAME/input"
NUMPREDICTIONS=40
NUMSAMPLES=5
DB_FILE="$BASEDIR/abag-benchmark-set/data/db/lightDb.txt"
MODELDIR="$BASEDIR/abag-benchmark-set/data/models"
FILTERED_PDB_DIR="$BASEDIR/abag-benchmark-set/data/db/structures_filtered"

# Templates for paths
PDB_TEMPLATE="$MODELDIR/$RUNNAME/{ID}complex/model_seed_{PRED_NUM}_sample_{MODEL_NUM}_${RUNNAME}.pdb"
FEATURES_TEMPLATE="$MODELDIR/$RUNNAME/{ID}complex/full_confidences_seed_{PRED_NUM}_sample_{MODEL_NUM}_${RUNNAME}.json"
DATA_TEMPLATE="$MODELDIR/$RUNNAME/{ID}complex/full_confidences_seed_{PRED_NUM}_sample_{MODEL_NUM}_${RUNNAME}.json"
SAMPLE_ID_TEMPLATE="{ID}_seed_{PRED_NUM_REFERENCE}_sample_{MODEL_NUM_REFERENCE}_vs_seed_{PRED_NUM_QUERY}_sample_{MODEL_NUM_QUERY}_${RUNNAME}"

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

substitute_sample_template() {
    local template="$1"
    local id="$2"
    local model_num_ref="$3"
    local pred_num_ref="$4"
    local model_num_query="$5"
    local pred_num_query="$6"

    template="${template//\{ID\}/$id}"
    template="${template//\{MODEL_NUM_REFERENCE\}/$model_num_ref}"
    template="${template//\{PRED_NUM_REFERENCE\}/$pred_num_ref}"
    template="${template//\{MODEL_NUM_QUERY\}/$model_num_query}"
    template="${template//\{PRED_NUM_QUERY\}/$pred_num_query}"
    echo "$template"
}

# Read the database file, skipping the header
tail -n +2 "$DB_FILE" | while IFS=, read -r pdb_id chain_A chain_H chain_L; do

    # Set output CSV file for the current pdb_id
    output_file="$OUTDIR/${pdb_id}_samples_$RUNNAME.csv"

    # Initialize the output file
    > "$output_file"
    echo "sample_id,pdbid,Achain,Hchain,Lchain,ground_truth_pdb,query_pdb,query_af_features,query_af_data,model_num_query,prediction_query,reference_pdb,reference_af_features,reference_af_data,model_num_reference,prediction_reference,preset" >> "$output_file"

    # Process all combinations of models and samples
    for ((MODEL_NUM_REF = 1; MODEL_NUM_REF <= NUMSAMPLES; MODEL_NUM_REF++)); do
        for ((PRED_NUM_REF = 0; PRED_NUM_REF < NUMPREDICTIONS; PRED_NUM_REF++)); do
            for ((MODEL_NUM_QUERY = 1; MODEL_NUM_QUERY <= NUMSAMPLES; MODEL_NUM_QUERY++)); do
                for ((PRED_NUM_QUERY = 0; PRED_NUM_QUERY < NUMPREDICTIONS; PRED_NUM_QUERY++)); do
                    
                    # Generate dynamic values
                    SAMPLE_ID=$(substitute_sample_template "$SAMPLE_ID_TEMPLATE" "$pdb_id" "$MODEL_NUM_REF" "$PRED_NUM_REF" "$MODEL_NUM_QUERY" "$PRED_NUM_QUERY")
                    
                    QUERY_PDB=$(substitute_template "$PDB_TEMPLATE" "$pdb_id" "$MODEL_NUM_QUERY" "$PRED_NUM_QUERY")
                    QUERY_DATA=$(substitute_template "$DATA_TEMPLATE" "$pdb_id" "$MODEL_NUM_QUERY" "$PRED_NUM_QUERY")
                    QUERY_FEATURES=$(substitute_template "$FEATURES_TEMPLATE" "$pdb_id" "$MODEL_NUM_QUERY" "$PRED_NUM_QUERY")

                    REF_PDB=$(substitute_template "$PDB_TEMPLATE" "$pdb_id" "$MODEL_NUM_REF" "$PRED_NUM_REF")
                    REF_DATA=$(substitute_template "$DATA_TEMPLATE" "$pdb_id" "$MODEL_NUM_REF" "$PRED_NUM_REF")
                    REF_FEATURES=$(substitute_template "$FEATURES_TEMPLATE" "$pdb_id" "$MODEL_NUM_REF" "$PRED_NUM_REF")

                    GROUND_TRUTH_PDB="$FILTERED_PDB_DIR/${pdb_id}_filtered.pdb"

                    # Validate file existence
                    all_files_exist=true
                    for file in "$QUERY_PDB" "$REF_PDB" "$QUERY_DATA" "$REF_DATA" "$QUERY_FEATURES" "$REF_FEATURES"; do
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
                    echo "${SAMPLE_ID},${pdb_id},${chain_A},${chain_H},${chain_L},${GROUND_TRUTH_PDB},${QUERY_PDB},${QUERY_FEATURES},${QUERY_DATA},${MODEL_NUM_QUERY},${PRED_NUM_QUERY},${REF_PDB},${REF_FEATURES},${REF_DATA},${MODEL_NUM_REF},${PRED_NUM_REF},${RUNNAME}" >> "$output_file"
                done
            done
        done
    done
done

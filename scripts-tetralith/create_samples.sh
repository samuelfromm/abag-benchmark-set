#!/bin/bash

# Check if BASEDIR argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <RUNNAME>"
    exit 1
fi

# Set BASEDIR from the command-line argument
BASEDIR="/proj/elofssonlab/users/x_safro/git"

# Set the default values for the variables
RUNNAME="$1"
NUMPREDICTIONS=40
WEIGHT_VERSION="multimer_v3"
DB_FILE="$BASEDIR/abag-benchmark-set/data/db/lightDb.txt"
MODELDIR="$BASEDIR/abag-benchmark-set/data/models"

OUTPUT_FILE="./samples_$RUNNAME.csv"
FILTERED_PDB_DIR="$BASEDIR/abag-benchmark-set/data/db/structures_filtered"

# Templates for paths
PDB_TEMPLATE="$MODELDIR/$RUNNAME/{ID}complex/unrelaxed_model_{MODEL_NUM}_${WEIGHT_VERSION}_pred_{PRED_NUM}_${RUNNAME}.pdb"
FEATURES_TEMPLATE="$MODELDIR/$RUNNAME/{ID}complex/features.pkl"
DATA_TEMPLATE="$MODELDIR/$RUNNAME/{ID}complex/result_model_{MODEL_NUM}_${WEIGHT_VERSION}_pred_{PRED_NUM}_${RUNNAME}.pkl"
SAMPLE_ID_TEMPLATE="{ID}_model_{MODEL_NUM}_${WEIGHT_VERSION}_pred_{PRED_NUM}_${RUNNAME}"

# Initialize the output file
> "$OUTPUT_FILE"
echo "sample_id,pdbid,Achain,Hchain,Lchain,reference_pdb,query_pdb,query_af_features,query_af_data,model,prediction" >> "$OUTPUT_FILE"

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
    # if [[ "$pdb_id" != "7ox2" ]]; then
    #     continue
    # fi
    for MODEL_NUM in {1..5}; do
        for ((PRED_NUM = 0; PRED_NUM < NUMPREDICTIONS; PRED_NUM++)); do

            # Generate dynamic values
            SAMPLE_ID=$(substitute_template "$SAMPLE_ID_TEMPLATE" "$pdb_id" "$MODEL_NUM" "$PRED_NUM")
            QUERY_PDB=$(substitute_template "$PDB_TEMPLATE" "$pdb_id" "$MODEL_NUM" "$PRED_NUM")
            QUERY_DATA=$(substitute_template "$DATA_TEMPLATE" "$pdb_id" "$MODEL_NUM" "$PRED_NUM")
            QUERY_FEATURES=$(substitute_template "$FEATURES_TEMPLATE" "$pdb_id" "$MODEL_NUM" "$PRED_NUM")
            REFERENCE_PDB="$FILTERED_PDB_DIR/${pdb_id}_filtered.pdb"
            MODEL="model_${MODEL_NUM}_${WEIGHT_VERSION}"

            # Validate file existence and skip appending to OUTPUT_FILE if any file is missing
            all_files_exist=true
            for file in "$QUERY_PDB" "$REFERENCE_PDB" "$QUERY_DATA" "$QUERY_FEATURES"; do
                if [ ! -f "$file" ]; then
                    echo "Warning: File $file not found"
                    all_files_exist=false
                    break
                fi
            done

            # Skip processing if any file is missing
            if [ "$all_files_exist" = false ]; then
                continue
            fi

            # Append the processed data to the output file
            echo "${SAMPLE_ID},${pdb_id},${chain_A},${chain_H},${chain_L},${REFERENCE_PDB},${QUERY_PDB},${QUERY_FEATURES},${QUERY_DATA},${MODEL},${PRED_NUM}" >> "$OUTPUT_FILE"
        done
    done
done
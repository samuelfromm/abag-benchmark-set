#!/bin/bash

# Check if RUNNAME argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <RUNNAME>"
    exit 1
fi

# Set variables
RUNNAME="$1"
BASEDIR="/proj/berzelius-2021-29/users/x_safro/git/paper"
OUTDIR="$BASEDIR/abag-benchmark-set/data/scores_paired/$RUNNAME/input"
SOURCE_DIR="$BASEDIR/abag-benchmark-set/data/scores_ae/$RUNNAME/output"
NUMPREDICTIONS=40
NUMSAMPLES=5
NSAMPLES=400 # number of random symmetric-unique samples per pdb_id
DB_FILE="$BASEDIR/abag-benchmark-set/data/db/lightDb.txt"
MODELDIR="$BASEDIR/abag-benchmark-set/data/models"
FILTERED_PDB_DIR="$BASEDIR/abag-benchmark-set/data/db/structures_filtered"

# Templates
PDB_TEMPLATE="$MODELDIR/$RUNNAME/{ID}complex/model_seed_{PRED_NUM}_sample_{MODEL_NUM}_${RUNNAME}.pdb"
FEATURES_TEMPLATE="$MODELDIR/$RUNNAME/{ID}complex/full_confidences_seed_{PRED_NUM}_sample_{MODEL_NUM}_${RUNNAME}.json"
DATA_TEMPLATE="$MODELDIR/$RUNNAME/{ID}complex/full_confidences_seed_{PRED_NUM}_sample_{MODEL_NUM}_${RUNNAME}.json"
SAMPLE_ID_TEMPLATE="{ID}_seed_{PRED_NUM_REFERENCE}_sample_{MODEL_NUM_REFERENCE}_vs_seed_{PRED_NUM_QUERY}_sample_{MODEL_NUM_QUERY}_${RUNNAME}"
SAMPLE_ID_TEMPLATE_SINGLE="{ID}_seed_{PRED_NUM}_sample_{MODEL_NUM}_${RUNNAME}"
ABAG_PDB_TEMPLATE="$SOURCE_DIR/{ID}/{ID}_seed_{PRED_NUM}_sample_{MODEL_NUM}_${RUNNAME}/{ID}_seed_{PRED_NUM}_sample_{MODEL_NUM}_${RUNNAME}_query_abag_merged.cif"
AE_TEMPLATE="$SOURCE_DIR/{ID}/{ID}_seed_{PRED_NUM}_sample_{MODEL_NUM}_${RUNNAME}/{ID}_seed_{PRED_NUM}_sample_{MODEL_NUM}_${RUNNAME}_aligned_error.pth"

# Template substitution functions
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

# Generate one random non-self symmetric-invariant pair
generate_random_pair() {
    local m1 p1 m2 p2
    while true; do
        m1=$(( RANDOM % NUMSAMPLES + 1 ))
        p1=$(( RANDOM % NUMPREDICTIONS ))
        m2=$(( RANDOM % NUMSAMPLES + 1 ))
        p2=$(( RANDOM % NUMPREDICTIONS ))

        # Exclude trivial identical pairs
        if [[ $m1 -eq $m2 && $p1 -eq $p2 ]]; then
            continue
        fi

        echo "$m1 $p1 $m2 $p2"
        return
    done
}

# Symmetry normalization key
make_key() {
    local a b c d
    a="$1"; b="$2"; c="$3"; d="$4"

    if [[ $a -lt $c ]]; then
        echo "$a,$b,$c,$d"
    elif [[ $a -gt $c ]]; then
        echo "$c,$d,$a,$b"
    else
        if [[ $b -le $d ]]; then
            echo "$a,$b,$c,$d"
        else
            echo "$c,$d,$a,$b"
        fi
    fi
}

# Read the database file

MAX_TRIES=5000


mkdir -p "$OUTDIR"

#tail -n +2 "$DB_FILE" | while IFS=, read -r pdb_id chain_A chain_H chain_L; do
while IFS=, read -r pdb_id chain_A chain_H chain_L; do

    echo "Processing PDB ID: $pdb_id"
    output_file="$OUTDIR/${pdb_id}_samples_$RUNNAME.csv"

    > "$output_file"
    echo "sample_id,sample_id_ref,sample_id_query,pdbid,Achain,Hchain,Lchain,ground_truth_pdb,query_pdb,query_af_features,query_af_data,model_num_query,prediction_query,reference_pdb,reference_af_features,reference_af_data,model_num_reference,prediction_reference,preset,query_pdb_abag,query_ae,reference_pdb_abag,reference_ae" >> "$output_file"

    unset seen
    declare -A seen=()
    count=0
    tries=0
    while [[ $count -lt $NSAMPLES && $tries -lt $MAX_TRIES ]]; do

        ((tries++))
        if [[ $tries -ge $MAX_TRIES ]]; then
            echo "Reached maximum number of tries ($MAX_TRIES) for PDB ID: $pdb_id"
            break
        fi
        read mref pref mquery pquery < <(generate_random_pair)

        key=$(make_key "$mref" "$pref" "$mquery" "$pquery")
        if [[ -n "${seen[$key]}" ]]; then
            continue
        fi
        seen[$key]=1

        SAMPLE_ID=$(substitute_sample_template "$SAMPLE_ID_TEMPLATE" "$pdb_id" "$mref" "$pref" "$mquery" "$pquery")
        SAMPLE_ID_REF=$(substitute_template "$SAMPLE_ID_TEMPLATE_SINGLE" "$pdb_id" "$mref" "$pref")
        SAMPLE_ID_QUERY=$(substitute_template "$SAMPLE_ID_TEMPLATE_SINGLE" "$pdb_id" "$mquery" "$pquery")

        QUERY_PDB=$(substitute_template "$PDB_TEMPLATE" "$pdb_id" "$mquery" "$pquery")
        QUERY_DATA=$(substitute_template "$DATA_TEMPLATE" "$pdb_id" "$mquery" "$pquery")
        QUERY_FEATURES=$(substitute_template "$FEATURES_TEMPLATE" "$pdb_id" "$mquery" "$pquery")
        QUERY_PDB_ABAG=$(substitute_template "$ABAG_PDB_TEMPLATE" "$pdb_id" "$mquery" "$pquery")
        QUERY_AE=$(substitute_template "$AE_TEMPLATE" "$pdb_id" "$mquery" "$pquery")

        REF_PDB=$(substitute_template "$PDB_TEMPLATE" "$pdb_id" "$mref" "$pref")
        REF_DATA=$(substitute_template "$DATA_TEMPLATE" "$pdb_id" "$mref" "$pref")
        REF_FEATURES=$(substitute_template "$FEATURES_TEMPLATE" "$pdb_id" "$mref" "$pref")
        REF_PDB_ABAG=$(substitute_template "$ABAG_PDB_TEMPLATE" "$pdb_id" "$mref" "$pref")
        REF_AE=$(substitute_template "$AE_TEMPLATE" "$pdb_id" "$mref" "$pref")

        GROUND_TRUTH_PDB="$FILTERED_PDB_DIR/${pdb_id}_filtered.pdb"

        all_files_exist=true
        for file in "$QUERY_PDB" "$REF_PDB" "$QUERY_DATA" "$REF_DATA" "$QUERY_FEATURES" "$REF_FEATURES" "$QUERY_PDB_ABAG" "$REF_PDB_ABAG" "$QUERY_AE" "$REF_AE" "$GROUND_TRUTH_PDB"; do
            if [ ! -f "$file" ]; then
                echo "Warning: File $file not found for $pdb_id"
                all_files_exist=false
                ((tries++))
                tries=$MAX_TRIES
                break
            fi
        done
        if ! $all_files_exist; then
            continue
        fi

        echo "${SAMPLE_ID},${SAMPLE_ID_REF},${SAMPLE_ID_QUERY},${pdb_id},${chain_A},${chain_H},${chain_L},${GROUND_TRUTH_PDB},${QUERY_PDB},${QUERY_FEATURES},${QUERY_DATA},${mquery},${pquery},${REF_PDB},${REF_FEATURES},${REF_DATA},${mref},${pref},${RUNNAME},${QUERY_PDB_ABAG},${QUERY_AE},${REF_PDB_ABAG},${REF_AE}" >> "$output_file"

        ((count++))
    done

done < <(tail -n +2 "$DB_FILE")
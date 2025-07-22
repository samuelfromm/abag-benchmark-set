#!/bin/bash

# Check for required arguments
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <ID_file.csv> <pdb_directory> <output_directory> [--no-header]"
    exit 1
fi

# Assign arguments to variables
ID_FILE="$1"
PDB_DIR="$2"
OUTPUT_DIR="$3"
NO_HEADER=false

# Check for the optional --no-header flag
if [ "$#" -eq 4 ] && [ "$4" == "--no-header" ]; then
    NO_HEADER=true
fi

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Define the start line (skip header if NO_HEADER is false)
START_LINE=1
if ! $NO_HEADER; then
    START_LINE=2
fi

# Read the CSV file line by line, handling each PDB ID and chain list
tail -n +$START_LINE "$ID_FILE" | while IFS=, read -r pdb_id chain_A chain_H chain_L; do
    # Construct the path to the input PDB file
    pdb_file="${PDB_DIR}/${pdb_id}.pdb"

    # Check if the PDB file exists
    if [ ! -f "$pdb_file" ]; then
        echo "Warning: PDB file for ID $pdb_id not found in $PDB_DIR"
        continue
    fi

    # Collect chains to keep (ignoring "NA" values) and handle multiple chains in each slot
    chains=()
    [[ "$chain_A" != "NA" ]] && chains+=($(echo "$chain_A" | tr '|' ' '))
    [[ "$chain_H" != "NA" ]] && chains+=($(echo "$chain_H" | tr '|' ' '))
    [[ "$chain_L" != "NA" ]] && chains+=($(echo "$chain_L" | tr '|' ' '))

    # If no valid chains are specified, skip this PDB file
    if [ "${#chains[@]}" -eq 0 ]; then
        echo "Warning: No valid chains specified for PDB ID $pdb_id"
        continue
    fi

    # Join chains with spaces for passing to the Python script
    chains_to_keep="${chains[@]}"

    # Define the output file path
    output_pdb_file="${OUTPUT_DIR}/${pdb_id}_filtered.pdb"

    # Run the Python script to filter the chains
    python3 filter_pdb_chains.py "$pdb_file" "$output_pdb_file" $chains_to_keep

    # Check if the Python script executed successfully
    if [ $? -eq 0 ]; then
        echo "Filtered PDB for $pdb_id saved to $output_pdb_file"
    else
        echo "Error: Failed to filter PDB for ID $pdb_id"
    fi
done
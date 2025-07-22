#!/bin/bash -x

# This script checks if the features.pkl file exists for each ID in the provided IDs file.
# If the file exists, the ID is appended to IDs_completed.csv.
# If the file does not exist, the ID is appended to IDs_missing.csv.

# Check if IDFILE and MSADIR are provided as command-line arguments
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <IDFILE> <MSADIR>"
    exit 1
fi

IDFILE=$1    # First argument is the path to the IDs file
MSADIR=$2    # Second argument is the directory containing MSA files

# Output files in the current directory
MISSINGFILE="./IDs_missing.csv"
COMPLETEDFILE="./IDs_completed.csv"

# Create empty output files
> "$MISSINGFILE"
> "$COMPLETEDFILE"

# Loop through each ID in the ID file
while IFS= read -r ID; do
    # Check if features.pkl exists in MSADIR/$ID
    if [ -f "$MSADIR/$ID/features.pkl" ]; then
        # Append ID to completed file
        echo "$ID" >> "$COMPLETEDFILE"
    else
        # Append ID to missing file
        echo "$ID" >> "$MISSINGFILE"
    fi
done < "$IDFILE"
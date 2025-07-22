#!/bin/bash -x

# This script combines all scores.csv files from directories listed in the provided IDs file
# into a single output file in the specified output directory.

# Check if IDFILE, SCOREDIR, and OUTDIR are provided as command-line arguments
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    echo "Usage: $0 <IDFILE> <SCOREDIR> <RUNNAME>"
    exit 1
fi

IDFILE=$1    # First argument is the path to the IDs file
SCOREDIR=$2  # Second argument is the directory containing scores files
RUNNAME=$3    # Third argument is the output directory for the combined file

# Output file
COMBINEDFILE="$SCOREDIR/$RUNNAME/scores_$RUNNAME.csv"

# Create or clear the output file and add a header
> "$COMBINEDFILE"

# Loop through each ID in the ID file
while IFS= read -r ID; do
    SCORE_FILE="$SCOREDIR/$RUNNAME/output/$ID/scores.csv"
    
    # Check if scores.csv exists for the current ID
    if [ -f "$SCORE_FILE" ]; then
        # Append the content of scores.csv to the combined file, excluding the header if it's not the first file
        if [ ! -s "$COMBINEDFILE" ]; then
            # First file: copy with header
            cat "$SCORE_FILE" >> "$COMBINEDFILE"
        else
            # Subsequent files: skip header
            tail -n +2 "$SCORE_FILE" >> "$COMBINEDFILE"
        fi
    fi
done < "$IDFILE"

echo "Scores combined into $COMBINEDFILE"
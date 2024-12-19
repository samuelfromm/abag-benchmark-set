#!/bin/bash -x

# This script checks if the features.pkl file exists for each ID in the provided IDs file.
# If the file exists, the ID is appended to IDs_completed.csv.
# If the file does not exist, the ID is appended to IDs_missing.csv.

# Check if IDFILE and MSADIR are provided as command-line arguments
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
    echo "Usage: $0 <IDFILE> <RUNNAME> <OUTPUTDIR> <NUMPREDICTIONS>"
    exit 1
fi

IDFILE=$1    # First argument is the path to the IDs file
RUNNAME=$2    # Second argument is the name of the run
OUTPUTDIR=$3    # Third argument is the output directory
NUMPREDICTIONS=$4    # Fourth argument is the number of predictions per model

# Output files in the current directory
MISSINGFILE="./IDs_${RUNNAME}_missing.csv"
COMPLETEDFILE="./IDs_${RUNNAME}_completed.csv"
ALLMISSINGFILE="./IDs_${RUNNAME}_all_missing.csv"

# Create empty output files
> "$MISSINGFILE"
> "$COMPLETEDFILE"
> "$ALLMISSINGFILE"

# Loop through each ID in the ID file
while IFS= read -r ID; do
    all_files_present=true
    # Check if all the files of type
    # OUTPUTDIR/RUNNAME/ID/result_model_{i}_multimer_v3_pred_{k}_{RUNNAME}.pkl
    # are present where i is 0 to 5 and k is 0 to NUMPREDICTIONS-1
    for i in {1..5}; do
        for ((k=0; k<NUMPREDICTIONS; k++)); do
            filepath="${OUTPUTDIR}/${RUNNAME}/${ID}/result_model_${i}_multimer_v3_pred_${k}_${RUNNAME}.pkl"
            if [ ! -f "$filepath" ]; then
                echo "Missing file: $filepath"
                all_files_present=false
                #break 2
                echo "${ID}/result_model_${i}_multimer_v3_pred_${k}_${RUNNAME}.pkl" >> "$ALLMISSINGFILE"
            fi
        done
    done

    if $all_files_present; then
        echo "$ID" >> "$COMPLETEDFILE"
    else
        echo "$ID" >> "$MISSINGFILE"
    fi
done < "$IDFILE"
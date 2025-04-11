import os
import json
import csv

# Function to process files in the directory
def process_files(base_dir, id_file):
    # Read IDs from the CSV file
    with open(id_file, 'r') as f:
        ids = [line.strip() for line in f if line.strip()]
    
    # Iterate over the IDs
    for id_ in ids:
        #id = id_+"complex"
        id = id_
        id_dir = os.path.join(base_dir, id)
        if not os.path.exists(id_dir):
            print(f"Directory for ID {id_} not found. Skipping...")
            continue

        # Process files in the ID directory
        for i in range(1,41):  # i from 0 to 40
            for j in range(5):  # j from 0 to 5
                pdb_filename = f"unrelaxed_model_{i}_multimer_v3_pred_{j}_alphafold3_{i}.pdb"
                result_json_filename = f"result_model_{i}_multimer_v3_pred_{j}_alphafold3_{i}.json"
                confidence_json_filename = f"full_confidences_{i}_multimer_v3_pred_{j}_alphafold3_{i}.json"

                pdb_path = os.path.join(id_dir, pdb_filename)
                result_json_path = os.path.join(id_dir, result_json_filename)
                confidence_json_path = os.path.join(id_dir, confidence_json_filename)

                # Check if all required files exist
                if not (os.path.exists(pdb_path) and os.path.exists(result_json_path) and os.path.exists(confidence_json_path)):
                    print(f"Missing files for i={i}, j={j} in ID {id_}. Skipping this set...")
                    continue

                # Rename the PDB file
                new_pdb_filename = f"model_seed_{i-1}_sample_{j+1}_alphafold3.pdb"
                new_pdb_path = os.path.join(id_dir, new_pdb_filename)
                os.rename(pdb_path, new_pdb_path)
                #print(f"Renamed {pdb_filename} to {new_pdb_filename}")

                # Combine JSON files
                new_combined_json_filename = f"full_confidences_seed_{i-1}_sample_{j+1}_alphafold3.json"
                new_combined_json_path = os.path.join(id_dir, new_combined_json_filename)
                try:
                    with open(result_json_path, 'r') as result_file, open(confidence_json_path, 'r') as confidence_file:
                        result_data = json.load(result_file)
                        confidence_data = json.load(confidence_file)
                        combined_data = {**result_data, **confidence_data}

                    with open(new_combined_json_path, 'w') as combined_file:
                        json.dump(combined_data, combined_file, indent=4)
                    #print(f"Combined JSON files into {new_combined_json_filename}")

                    # Delete original JSON files
                    os.remove(result_json_path)
                    os.remove(confidence_json_path)
                    #print(f"Deleted original JSON files for i={i}, j={j}")

                except Exception as e:
                    print(f"Error processing JSON files for i={i}, j={j}: {e}")

if __name__ == "__main__":
    # Set base directory and ID file path
    base_dir = "/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/models/alphafold3"  # Replace with your base directory
    id_file = "/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/db/IDs_msas.csv" #"/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/db/IDs_working.csv"  # Replace with the path to your CSV file

    # Call the processing function
    process_files(base_dir, id_file)
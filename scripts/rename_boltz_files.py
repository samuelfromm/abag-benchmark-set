import os
import json
import numpy as np
import pandas as pd
import shutil

# Path setup
base_dir = "/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/models/boltz0.4.0"
ids_csv_path = "/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/db/IDs_msas.csv"  # Update this to the path of your IDs CSV file
outdir = "/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/models/boltz"

# Load IDs from CSV
ids_df = pd.read_csv(ids_csv_path, header=None)  # Assuming the CSV has no headers
ids = ids_df[0].tolist()  # List of IDs

# Iterate through each ID and seed
for seed in range(1, 41):  # Seeds from 1 to 40
    for protein_id in ids:
        input_dir = os.path.join(base_dir, f"{protein_id}/output_{seed}/boltz_results_{protein_id}/predictions/{protein_id}")
        output_dir = os.path.join(outdir, f"{protein_id}/")  # Consolidated directory

        if not os.path.exists(input_dir):
            print(f"Directory {input_dir} does not exist. Skipping...")
            continue

        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Process files in the directory
        for i in range(5):  # Assuming model indices are from 0 to 4
            pdb_file = f"{protein_id}_model_{i}.pdb"
            json_file = f"confidence_{protein_id}_model_{i}.json"
            pae_file = f"pae_{protein_id}_model_{i}.npz"
            plddt_file = f"plddt_{protein_id}_model_{i}.npz"

            # Check if all required files exist
            if not all(os.path.exists(os.path.join(input_dir, f)) for f in [pdb_file, json_file, pae_file, plddt_file]):
                print(f"Missing files for {protein_id}, model {i}, seed {seed}. Skipping...")
                continue

            new_pdb_name = f"model_{i+1}_seed_{seed-1}_boltz.pdb"
            shutil.copy(
                os.path.join(input_dir, pdb_file),
                os.path.join(output_dir, new_pdb_name),
            )

            # Create new JSON file
            with open(os.path.join(input_dir, json_file), "r") as f:
                confidence_data = json.load(f)

            # Add new keys from NPZ files
            pae_data = np.load(os.path.join(input_dir, pae_file))["pae"].tolist()
            plddt_data = np.load(os.path.join(input_dir, plddt_file))["plddt"].tolist()

            confidence_data["pae"] = pae_data
            confidence_data["plddt"] = plddt_data

            new_json_name = f"full_confidence_model_{i+1}_seed_{seed-1}_boltz.json"
            with open(os.path.join(output_dir, new_json_name), "w") as f:
                json.dump(confidence_data, f, indent=4)


print("Processing complete.")
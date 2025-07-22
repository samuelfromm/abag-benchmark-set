import os
import json
import numpy as np
import pandas as pd
import torch
import logging
import shutil
from Bio.PDB import MMCIFParser, PDBIO

def convert_cif_to_pdb(input_cif_path, output_pdb_path):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("model", input_cif_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb_path)
    print(f"Converted: {input_cif_path} -> {output_pdb_path}")


logging.basicConfig(level=logging.INFO)

# Path setup
base_dir = "/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/models/chai1_raw"
ids_csv_path = "/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/db/IDs_msas.csv"  # Update this to the path of your IDs CSV file
outdir = "/proj/elofssonlab/users/x_safro/git/abag-benchmark-set/data/models/chai1"

# Load IDs from CSV
ids_df = pd.read_csv(ids_csv_path, header=None)  # Assuming the CSV has no headers
ids = ids_df[0].tolist()  # List of IDs

# # Keep only IDs before "8jn4complex"
# if "8jn4complex" in ids:
#     start_index = ids.index("8jn4complex")
#     ids = ids[:start_index]
# else:
#     ids = []  # Or raise an error if the element is expected

# Iterate through each ID and seed
for seed in range(1, 41):  # Seeds from 1 to 40
    logging.info(f"Working on seed {seed}")
    for protein_id in ids:
        logging.info(f"Working on id {protein_id}")
        input_dir = os.path.join(base_dir, f"seed{seed}/{protein_id}/")
        output_dir = os.path.join(outdir, f"{protein_id}/")  # Consolidated directory

        if not os.path.exists(input_dir):
            logging.warning(f"Directory {input_dir} does not exist. Skipping...")
            continue

        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Process files in the directory
        for i in range(5):  # Assuming model indices are from 0 to 4
            pdb_file = f"pred.model_idx_{i}.cif"
            pae_file = f"pae.pt"
            pde_file = f"pde.pt"
            plddt_file = f"plddt.pt"
            scores_file = f"scores.model_idx_{i}.npz"


            missing_files = False
            # Check if all required files exist
            for f in [pdb_file, pae_file, pde_file, plddt_file,scores_file]:
                if not os.path.exists(os.path.join(input_dir, f)):
                    logging.warning(f"Missing file: {f}")
                    logging.warning(f"Missing files for {protein_id}, model {i}, seed {seed}. Skipping...")
                    missing_files = True

            if missing_files:
                continue
            else:
                new_cif_name = f"model_{i+1}_seed_{seed-1}_chai1.cif"
                shutil.copy(
                    os.path.join(input_dir, pdb_file),
                    os.path.join(output_dir, new_cif_name),
                )

                new_pdb_file = f"model_{i+1}_seed_{seed-1}_chai1.pdb"
                convert_cif_to_pdb(os.path.join(input_dir, pdb_file),os.path.join(output_dir, new_pdb_file))

                # Add new keys from NPZ files
                data = np.load(os.path.join(input_dir, scores_file))

                # #KEYS
                # aggregate_score
                # ptm
                # iptm
                # per_chain_ptm
                # per_chain_pair_iptm
                # has_inter_chain_clashes
                # chain_chain_clashes

                plddt_data = torch.load(os.path.join(input_dir, plddt_file),weights_only=True)[i].numpy().tolist()
                pde_data = torch.load(os.path.join(input_dir, pde_file),weights_only=True)[i].numpy().tolist()
                pae_data = torch.load(os.path.join(input_dir, pae_file),weights_only=True)[i].numpy().tolist()

                confidence_data = {}
                confidence_data["pae"] = pae_data
                confidence_data["plddt"] = plddt_data
                confidence_data["pde"] = plddt_data
                for key in data:
                    confidence_data[key] = data[key][0].tolist()


                new_json_name = f"full_confidence_model_{i+1}_seed_{seed-1}_chai1.json"
                with open(os.path.join(output_dir, new_json_name), "w") as f:
                    json.dump(confidence_data, f, indent=4)
        logging.info(f"Finished id {protein_id}")


print("Processing complete.")
import argparse
import torch
import pandas as pd
import aligned_error.aligned_error_ops
import aligned_error.aligned_error_utils
import aligned_error.aligned_error_confidence
import pickle
import os
import numpy as np
import json
import aligned_error.confidence.confidence_tools

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Run aligned error confidence analysis.")
    parser.add_argument("--sample_id", required=True, help="Sample ID for processing.")
    parser.add_argument("--af_data_query", required=True, help="Path to the AF data pickle file for the query.")
    parser.add_argument("--af_features_query", required=True, help="Path to the AF features pickle file for the query.")
    parser.add_argument("--af_data_reference", required=True, help="Path to the AF data pickle file for the reference.")
    parser.add_argument("--af_features_reference", required=True, help="Path to the AF features pickle file for the reference.")
    parser.add_argument("--output_csv", required=True, help="Path to save the output CSV file.")
    return parser.parse_args()


def load_data(file_path):
    """Load data from a pickle file."""
    _, file_extension = os.path.splitext(file_path)
    
    if file_extension == ".pkl":
        # Load data from the pickle file
        with open(file_path, "rb") as f:
            data = pickle.load(f)
    elif file_extension == ".json":
        # Load data from the JSON file
        with open(file_path, "r") as f:
            data = json.load(f)
    else:
        print(f"Unsupported file format: {file_extension}. Please use a .pkl or .json file.")
        return
    
    # Print available keys in the data
    if not isinstance(data, dict):
        print("The data is not a dictionary.")
        return

    return data


def get_predicted_aligned_error(af_data, af_features):

    af_data = load_data(af_data)

    if "predicted_aligned_error" in af_data.keys():
        predicted_aligned_error_tensor = af_data['predicted_aligned_error']
    elif "pae" in af_data.keys():
        predicted_aligned_error_tensor = af_data["pae"]
    else:
        raise KeyError("Could not find the predicted aligned error data (keywords 'pae' or 'predicted_aligned_error').")

    # Convert predicted_aligned_error_tensor to a torch tensor
    if isinstance(predicted_aligned_error_tensor, list):
        predicted_aligned_error_tensor = torch.tensor(np.array(predicted_aligned_error_tensor), dtype=torch.float32)
    elif isinstance(predicted_aligned_error_tensor, np.ndarray):
        predicted_aligned_error_tensor = torch.tensor(predicted_aligned_error_tensor, dtype=torch.float32)
    elif not isinstance(predicted_aligned_error_tensor, torch.Tensor):
        raise TypeError("The 'pae' or 'predicted_aligned_error' value must be a list, NumPy array, or a Torch tensor.")

    features_af_data = load_data(af_features)

    if "asym_id" in features_af_data.keys():
        asym_id_from_features = features_af_data['asym_id']
    elif "token_chain_ids" in features_af_data.keys():
        asym_id_from_features = features_af_data["token_chain_ids"]
        # asym id is alphabetical, convert to numerical

        asym_id_mapping = {val: idx+1 for idx, val in enumerate(sorted(set(asym_id_from_features), key=asym_id_from_features.index))}
        asym_id_from_features = [asym_id_mapping[val] for val in asym_id_from_features]
    else:
        raise KeyError("Could not find the asymmetric id data (keywords 'asym_id' or 'token_chain_ids').")

    # Convert asym_id_from_features to a torch tensor
    if isinstance(asym_id_from_features, list):
        asym_id_from_features = torch.tensor(np.array(asym_id_from_features), dtype=torch.float32)
    elif isinstance(asym_id_from_features, np.ndarray):
        asym_id_from_features = torch.tensor(asym_id_from_features, dtype=torch.float32)
    elif not isinstance(asym_id_from_features, torch.Tensor):
        raise TypeError("The 'asym_id' or 'token_chain_ids' value must be a list, NumPy array, or a Torch tensor.")
    return predicted_aligned_error_tensor, asym_id_from_features

def main():
    args = parse_arguments()

    pae_query, asym_id_query = get_predicted_aligned_error(
        af_data=args.af_data_query,
        af_features=args.af_features_query,
    )

    pae_reference, asym_id_reference = get_predicted_aligned_error(
        af_data=args.af_data_reference,
        af_features=args.af_features_reference,
    )   

    assert torch.equal(asym_id_query, asym_id_reference), "Asym IDs from query and reference do not match."

    interface_contact_mask = aligned_error.aligned_error_ops.compute_interface_mask(
        asym_id=asym_id_query
    )

    mpae_interface = (torch.abs(pae_query - pae_reference))[interface_contact_mask == 1].mean()
    mpae = torch.abs(pae_query - pae_reference).mean()



    precision = 2
    # Prepare the output data
    output_data = {
        "sample_id": args.sample_id,
        "mpae_interface": round(mpae_interface.item(), precision),
        "mpae": round(mpae.item(), precision),
    }

    # Save the output data to a CSV file
    output_df = pd.DataFrame([output_data])
    output_df.to_csv(args.output_csv, index=False)


if __name__ == "__main__":
    main()

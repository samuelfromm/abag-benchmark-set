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

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Run predicted aligned error confidence analysis.")
    parser.add_argument("--sample_id", required=True, help="Sample ID for processing.")
    parser.add_argument("--af_data", required=True, help="Path to the AF data pickle file.")
    parser.add_argument("--af_features", required=True, help="Path to the AF features pickle file.")
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


def main():
    args = parse_arguments()
    precision = 2
    output_data = {
        "sample_id": args.sample_id,
    }

    af_data = load_data(args.af_data)

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

    features_af_data = load_data(args.af_features)

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

    # Calculate aeTM
    paetm = aligned_error.aligned_error_confidence.compute_aetm(
        aligned_error=predicted_aligned_error_tensor
    )

    # Calculate aeiTM
    paeitm = aligned_error.aligned_error_confidence.compute_aetm(
        aligned_error=predicted_aligned_error_tensor,
        asym_id=asym_id_from_features,
        interface=True,
    )

    pae_ranking_confidence = (
        aligned_error.aligned_error_confidence.calculate_ranking_confidence(
            ptm=paetm, iptm=paeitm
        )
    )

    # Calculate aeiTM per chain
    per_chain_aeitm = aligned_error.aligned_error_confidence.compute_per_chain_aetm(
        aligned_error=predicted_aligned_error_tensor,
        asym_id=asym_id_from_features,
    )
    per_chain_aeitm = [val.item() for val in per_chain_aeitm.values()]

    mean_per_chain_aeitm = sum(per_chain_aeitm) / len(per_chain_aeitm)
    max_per_chain_aeitm = max(per_chain_aeitm)
    min_per_chain_aeitm = min(per_chain_aeitm)

    # Prepare the output data
    output_data.update(
        {
            "paetm": round(paetm.item(), precision),
            "paeitm": round(paeitm.item(), precision),
            "pae_ranking_confidence": round(pae_ranking_confidence.item(), precision),
            "max_per_chain_paeitm": round(max_per_chain_aeitm, precision),
            "min_per_chain_paeitm": round(min_per_chain_aeitm, precision),
            "mean_per_chain_paeitm": round(mean_per_chain_aeitm, precision),
        }
    )

    # Save the output data to a CSV file
    output_df = pd.DataFrame([output_data])
    output_df.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    main()

import argparse
import torch
import pandas as pd
import pickle
import os
import json

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Process AlphaFold data and extract statistics.")
    parser.add_argument("--sample_id", required=True, help="Sample ID for processing.")
    parser.add_argument("--af_data", required=True, help="Path to the AlphaFold data pickle file.")
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



def process_af_data(af_data, precision=2):
    """Extract relevant metrics from AlphaFold data, handling tensors and floats/ints."""
    
    def safe_extract(value, precision):
        """Extracts the value from a tensor or returns it directly if not a tensor."""
        if isinstance(value, torch.Tensor):
            return round(value.item(), precision)
        return round(value, precision)

    ptm = af_data["ptm"]
    iptm = af_data["iptm"]
    
    if "ranking_score" in af_data:
        ranking_confidence = af_data["ranking_score"]
    else:
        ranking_confidence = af_data["ranking_confidence"]
    
    if "num_recycles" in af_data:
        num_recycles = af_data["num_recycles"]
    else:
        num_recycles = 0

    return {
        "ptm": safe_extract(ptm, precision),
        "iptm": safe_extract(iptm, precision),
        "ranking_confidence": safe_extract(ranking_confidence, precision),
        "num_recycles": safe_extract(num_recycles, 0)  
    }

def save_to_csv(output_data, output_csv_path):
    """Save output data to a CSV file."""
    output_df = pd.DataFrame([output_data])
    output_df.to_csv(output_csv_path, index=False)

def main():
    args = parse_arguments()
    af_data = load_data(args.af_data)
    extracted_metrics = process_af_data(af_data)

    output_data = {
        "sample_id": args.sample_id,
        **extracted_metrics,
    }

    save_to_csv(output_data, args.output_csv)

if __name__ == "__main__":
    main()
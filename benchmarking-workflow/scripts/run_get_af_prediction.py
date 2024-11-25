import argparse
import torch
import pandas as pd
import pickle

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Process AlphaFold data and extract statistics.")
    parser.add_argument("--sample_id", required=True, help="Sample ID for processing.")
    parser.add_argument("--af_data", required=True, help="Path to the AlphaFold data pickle file.")
    parser.add_argument("--output_csv", required=True, help="Path to save the output CSV file.")
    return parser.parse_args()

def load_data_from_pkl(pkl_path):
    """Load data from a pickle file."""
    with open(pkl_path, "rb") as file:
        return pickle.load(file)

def process_af_data(af_data, precision=2):
    """Extract relevant metrics from AlphaFold data."""
    ptm = torch.from_numpy(af_data["ptm"])
    iptm = torch.from_numpy(af_data["iptm"])
    ranking_confidence = torch.tensor(af_data["ranking_confidence"])  # Not a numpy array, but a float
    num_recycles = torch.tensor(af_data["num_recycles"])

    return {
        "ptm": round(ptm.item(), precision),
        "iptm": round(iptm.item(), precision),
        "ranking_confidence": round(ranking_confidence.item(), precision),
        "num_recycles": num_recycles.item(),
    }

def save_to_csv(output_data, output_csv_path):
    """Save output data to a CSV file."""
    output_df = pd.DataFrame([output_data])
    output_df.to_csv(output_csv_path, index=False)

def main():
    args = parse_arguments()
    af_data = load_data_from_pkl(args.af_data)
    extracted_metrics = process_af_data(af_data)

    output_data = {
        "sample_id": args.sample_id,
        **extracted_metrics,
    }

    save_to_csv(output_data, args.output_csv)

if __name__ == "__main__":
    main()
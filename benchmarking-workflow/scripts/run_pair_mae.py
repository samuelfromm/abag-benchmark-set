import argparse
import torch
import pandas as pd
import aligned_error.aligned_error_ops
import aligned_error.aligned_error_utils
import aligned_error.aligned_error_confidence

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Run aligned error confidence analysis.")
    parser.add_argument("--sample_id", required=True, help="Sample ID for processing.")
    parser.add_argument("--query_aligned_error", required=True, help="Path to the query aligned error file.")
    parser.add_argument("--reference_aligned_error", required=True, help="Path to the reference aligned error file.")
    parser.add_argument("--output_csv", required=True, help="Path to save the output CSV file.")
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    query_aligned_error = torch.load(args.query_aligned_error)
    query_ae = query_aligned_error["aligned_error"]
    query_asym_id = query_aligned_error["asym_id"]


    reference_aligned_error = torch.load(args.reference_aligned_error)
    reference_ae = reference_aligned_error["aligned_error"]
    reference_asym_id = reference_aligned_error["asym_id"]

    assert torch.equal(reference_asym_id, query_asym_id), "Asym IDs from query and reference do not match."


    interface_contact_mask = aligned_error.aligned_error_ops.compute_interface_mask(
        asym_id=query_asym_id
    )


    mae_interface = (torch.abs(query_ae - reference_ae))[interface_contact_mask == 1].mean()
    mae = torch.abs(query_ae - reference_ae).mean()



    precision = 2
    # Prepare the output data
    output_data = {
        "sample_id": args.sample_id,
        "mae_interface": round(mae_interface.item(), precision),
        "mae": round(mae.item(), precision),
    }

    # Save the output data to a CSV file
    output_df = pd.DataFrame([output_data])
    output_df.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    main()

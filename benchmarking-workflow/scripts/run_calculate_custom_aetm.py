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
    parser.add_argument("--aligned_error_data", required=True, help="Path to the input aligned error data file.")
    parser.add_argument("--reference_pdb", required=True, help="Path to the reference PDB file.")
    parser.add_argument("--query_pdb", required=True, help="Path to the query PDB file.")
    parser.add_argument("--output_csv", required=True, help="Path to save the output CSV file.")
    return parser.parse_args()

def main():
    args = parse_arguments()

    aligned_error_data = torch.load(args.aligned_error_data)
    aligned_error_tensor = aligned_error_data["aligned_error"]
    asym_id = aligned_error_data["asym_id"]

    assert aligned_error_data["reference_pdb"] == args.reference_pdb
    assert aligned_error_data["query_pdb"] == args.query_pdb

    # Calculate aeTM
    aetm = aligned_error.aligned_error_confidence.compute_custom_aetm(
        aligned_error=aligned_error_tensor
    )

    # Calculate aeiTM
    aeitm = aligned_error.aligned_error_confidence.compute_custom_aetm(
        aligned_error=aligned_error_tensor,
        asym_id=asym_id,
        interface=True,
    )

    ae_ranking_confidence = (
        aligned_error.aligned_error_confidence.calculate_ranking_confidence(
            ptm=aetm, iptm=aeitm
        )
    )

    # Calculate aeiTM per chain
    per_chain_aeitm = aligned_error.aligned_error_confidence.compute_per_chain_custom_aetm(
        aligned_error=aligned_error_tensor,
        asym_id=asym_id,
    )
    per_chain_aeitm = [val.item() for val in per_chain_aeitm.values()]

    mean_per_chain_aeitm = sum(per_chain_aeitm) / len(per_chain_aeitm)
    max_per_chain_aeitm = max(per_chain_aeitm)
    min_per_chain_aeitm = min(per_chain_aeitm)

    precision = 2
    # Prepare the output data
    output_data = {
        "sample_id": args.sample_id,
        "custom_aetm": round(aetm.item(), precision),
        "custom_aeitm": round(aeitm.item(), precision),
        "custom_ae_ranking_confidence": round(ae_ranking_confidence.item(), precision),
        "max_per_chain_custom_aeitm": round(max_per_chain_aeitm, precision),
        "min_per_chain_custom_aeitm": round(min_per_chain_aeitm, precision),
        "mean_per_chain_custom_aeitm": round(mean_per_chain_aeitm, precision),
    }

    # Save the output data to a CSV file
    output_df = pd.DataFrame([output_data])
    output_df.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    main()

import argparse
import aligned_error.aligned_error_ops
import aligned_error.aligned_error_utils
import torch

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Run aligned error analysis with merged chains.")
    parser.add_argument("--sample_id", required=True, help="Sample ID for processing.")
    parser.add_argument("--query_pdb", required=True, help="Path to the query PDB file.")
    parser.add_argument("--reference_pdb", required=True, help="Path to the reference PDB file.")
    parser.add_argument("--output_data", required=True, help="Path to save the output file.")
    return parser.parse_args()

def main():
    args = parse_arguments()

    # Read the reference and query PDB structures
    reference_protein = aligned_error.aligned_error_utils.protein_from_pdb(args.reference_pdb)
    query_protein = aligned_error.aligned_error_utils.protein_from_pdb(args.query_pdb)

    # Calculate aligned error between the reference and query structures
    aligned_error_value, asym_id = aligned_error.aligned_error_ops.calculate_aligned_error(
        reference_protein=reference_protein,
        query_protein=query_protein,
        return_asym_id=True,
    )

    # Prepare output data
    output_data = {
        "sample_id": args.sample_id,
        "reference_pdb": args.reference_pdb,
        "query_pdb": args.query_pdb,
        "aligned_error": aligned_error_value,
        "asym_id": asym_id,
    }

    # Save the output data using torch
    torch.save(output_data, args.output_data)

if __name__ == "__main__":
    main()
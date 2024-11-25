import argparse
import bioutils.align_cut_renumber
import pandas as pd


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Perform alignment, cut, and renumbering.")
    parser.add_argument("--alignment", required=True, help="Path to the alignment CSV file.")
    parser.add_argument("--reference_pdb", required=True, help="Path to the reference PDB file.")
    parser.add_argument("--query_pdb", required=True, help="Path to the query PDB file.")
    parser.add_argument("--sample_id", required=True, help="Sample ID to process.")
    parser.add_argument("--reference_cut", required=True, help="Path to the output reference cut PDB file.")
    parser.add_argument("--query_cut", required=True, help="Path to the output query cut PDB file.")
    parser.add_argument("--output", required=True, help="Path to the output CSV file.")
    return parser.parse_args()



def main():
    args = parse_arguments()

    # Load the alignment data
    alignment_df = pd.read_csv(args.alignment, index_col="sample_id")

    # Extract the alignment for the current sample
    sample_id = args.sample_id
    alignment = [
        alignment_df.loc[sample_id, "aln_reference"],
        alignment_df.loc[sample_id, "aln_query"],
    ]

    # Perform alignment, cut, and renumbering
    (
        _,
        _,
        fraction_aligned_residues_native,
        fraction_aligned_residues_model,
        min_fraction_aligned_chain_residues_native,
        min_fraction_aligned_chain_residues_model,
    ) = bioutils.align_cut_renumber.align_and_cut(
        path_to_native=args.reference_pdb,
        path_to_model=args.query_pdb,
        path_native_out=args.reference_cut,
        path_model_out=args.query_cut,
        alignment=alignment,
        rename_chains_rule="native",
        return_string=False,
    )

    # Prepare the output data
    output_data = {
        "sample_id": sample_id,
        "fraction_aligned_residues_reference": fraction_aligned_residues_native,
        "fraction_aligned_residues_query": fraction_aligned_residues_model,
        "min_fraction_aligned_chain_residues_reference": min_fraction_aligned_chain_residues_native,
        "min_fraction_aligned_chain_residues_query": min_fraction_aligned_chain_residues_model,
        "reference_cut": args.reference_cut,
        "query_cut": args.query_cut,
    }

    # Save the output data to a CSV file
    output_df = pd.DataFrame([output_data])
    output_df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()
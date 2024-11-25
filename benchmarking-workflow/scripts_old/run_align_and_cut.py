from snakemake.script import snakemake
import bioutils.align_cut_renumber
import pandas as pd

# Load the alignment data
alignment_df = pd.read_csv(snakemake.input.alignment, index_col="sample_id")

# Extract the alignment for the current sample
sample_id = snakemake.wildcards.sample_id
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
    path_to_native=snakemake.input.reference_pdb,
    path_to_model=snakemake.input.query_pdb,
    path_native_out=snakemake.output.reference_cut,
    path_model_out=snakemake.output.query_cut,
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
    "reference_cut": snakemake.output.reference_cut,
    "query_cut": snakemake.output.query_cut,
}

# Save the output data to a CSV file
output_df = pd.DataFrame([output_data])
output_df.to_csv(snakemake.output[0], index=False)

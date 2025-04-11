from snakemake.script import snakemake
import DockQ.DockQ
import pandas as pd

# Load the alignment data
alignment_df = pd.read_csv(snakemake.input.alignment, index_col="sample_id")

# Extract the alignment for the current sample
sample_id = snakemake.wildcards.sample_id
alignment_reference = alignment_df.loc[sample_id, "aln_reference_cut"]
alignment_query = alignment_df.loc[sample_id, "aln_query_cut"]

# Split the alignment data into chain IDs
native_chain_ids = [chain_id for chain_id in alignment_reference.split(":") if chain_id]
model_chain_ids = [chain_id for chain_id in alignment_query.split(":") if chain_id]

# Create a chain map between model and native interfaces
chain_map = dict(zip(model_chain_ids, native_chain_ids))

# Load the model and native structures
model = DockQ.DockQ.load_PDB(snakemake.input.query_cut)
native = DockQ.DockQ.load_PDB(snakemake.input.reference_cut)

# Run DockQ and get the output data
dockq_results, total_dockq_score = DockQ.DockQ.run_on_all_native_interfaces(
    model, native, chain_map=chain_map
)

# Initialize an empty dictionary to store statistics for each metric
metrics = ["DockQ", "LRMSD", "iRMSD", "fnat", "clashes"]
output_data = {"sample_id": sample_id}
precision = 2

# Loop through each metric and calculate statistics
for metric in metrics:
    # Extract metric values for all interfaces
    metric_vals = [interface_data[metric] for interface_data in dockq_results.values()]

    # Calculate statistics for the current metric
    average_metric = sum(metric_vals) / len(metric_vals)
    max_metric = max(metric_vals)
    min_metric = min(metric_vals)

    # Store the results in the output data with rounded values
    output_data[f"min_{metric.lower()}"] = round(min_metric, precision)
    output_data[f"max_{metric.lower()}"] = round(max_metric, precision)
    output_data[f"average_{metric.lower()}"] = round(average_metric, precision)


# Save the output data to a CSV file
output_df = pd.DataFrame([output_data])
output_df.to_csv(snakemake.output[0], index=False)


# ----------------------EXAMPLE-DOCKQ--------------------#
# (
#     {
#         "BC": {
#             "DockQ": 0.8244666441924884,
#             "F1": 0.9737827715355806,
#             "iRMSD": 1.4156530786424455,
#             "LRMSD": 1.9059374687955182,
#             "fnat": 0.9923664122137404,
#             "nat_correct": 130,
#             "nat_total": 131,
#             "fnonnat": 0.04411764705882353,
#             "nonnat_count": 6,
#             "model_total": 136,
#             "clashes": 1,
#             "len1": 221,
#             "len2": 215,
#             "class1": "receptor",
#             "class2": "ligand",
#             "is_het": False,
#             "chain1": "B",
#             "chain2": "C",
#             "chain_map": {"B": "B", "C": "C", "D": "D"},
#         },
#         "BD": {
#             "DockQ": 0.471769111353139,
#             "F1": 0.5,
#             "iRMSD": 2.5348365669333424,
#             "LRMSD": 8.399391122373226,
#             "fnat": 0.65,
#             "nat_correct": 13,
#             "nat_total": 20,
#             "fnonnat": 0.59375,
#             "nonnat_count": 19,
#             "model_total": 32,
#             "clashes": 2,
#             "len1": 221,
#             "len2": 197,
#             "class1": "receptor",
#             "class2": "ligand",
#             "is_het": False,
#             "chain1": "B",
#             "chain2": "D",
#             "chain_map": {"B": "B", "C": "C", "D": "D"},
#         },
#         "CD": {
#             "DockQ": 0.3336184813685321,
#             "F1": 0.3076923076923077,
#             "iRMSD": 2.480613029682055,
#             "LRMSD": 9.87845598600075,
#             "fnat": 0.3076923076923077,
#             "nat_correct": 4,
#             "nat_total": 13,
#             "fnonnat": 0.6923076923076923,
#             "nonnat_count": 9,
#             "model_total": 13,
#             "clashes": 1,
#             "len1": 215,
#             "len2": 197,
#             "class1": "receptor",
#             "class2": "ligand",
#             "is_het": False,
#             "chain1": "C",
#             "chain2": "D",
#             "chain_map": {"B": "B", "C": "C", "D": "D"},
#         },
#     },
#     1.6298542369141595,
# )

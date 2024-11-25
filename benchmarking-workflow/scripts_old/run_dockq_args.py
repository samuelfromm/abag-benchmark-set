import argparse
import DockQ.DockQ
import pandas as pd

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Process DockQ alignment data.")
parser.add_argument("--alignment", required=True, help="Path to the alignment CSV file")
parser.add_argument("--query_cut", required=True, help="Path to the query PDB file")
parser.add_argument("--reference_cut", required=True, help="Path to the reference PDB file")
parser.add_argument("--sample_id", required=True, help="Sample ID to process")
parser.add_argument("--output", required=True, help="Path to the output CSV file")

args = parser.parse_args()

# Load the alignment data
alignment_df = pd.read_csv(args.alignment, index_col="sample_id")

# Extract the alignment for the current sample
sample_id = args.sample_id
alignment_reference = alignment_df.loc[sample_id, "aln_reference_cut"]
alignment_query = alignment_df.loc[sample_id, "aln_query_cut"]

# Split the alignment data into chain IDs
native_chain_ids = [chain_id for chain_id in alignment_reference.split(":") if chain_id]
model_chain_ids = [chain_id for chain_id in alignment_query.split(":") if chain_id]

# Create a chain map between model and native interfaces
chain_map = dict(zip(model_chain_ids, native_chain_ids))

# Load the model and native structures
model = DockQ.DockQ.load_PDB(args.query_cut)
native = DockQ.DockQ.load_PDB(args.reference_cut)

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
output_df.to_csv(args.output, index=False)
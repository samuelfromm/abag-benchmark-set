from snakemake.script import snakemake
import DockQ.DockQ
import pandas as pd

def merge_chains(model, chains_to_merge, new_chain_name=None):
    for chain in chains_to_merge[1:]:
        for i, res in enumerate(model[chain]):
            res.id = (chain, res.id[1], res.id[2])
            model[chains_to_merge[0]].add(res)
        model.detach_child(chain)
    if new_chain_name:
        model[chains_to_merge[0]].id = "".join(new_chain_name)
    else:
        model[chains_to_merge[0]].id = "".join(chains_to_merge)
    return model



sample_id = snakemake.wildcards.sample_id

input_samples_df = pd.read_csv(snakemake.input[0], index_col="sample_id", na_values=[], keep_default_na=False)

Hchain = input_samples_df.at[sample_id, "Hchain"]
Achain = input_samples_df.at[sample_id, "Achain"]
Lchain = input_samples_df.at[sample_id, "Lchain"]


agchains = [part.strip() for part in Achain.split("|")]
abchains = [Hchain]
if Lchain != "NA":
    abchains.append(Lchain)

print("Ag chains: ", agchains)
print("Ab chains: ", abchains)

# Load the model and native structures
model = DockQ.DockQ.load_PDB(snakemake.input.query_cut)
native = DockQ.DockQ.load_PDB(snakemake.input.reference_cut)

# merge chains and rename chains

model = merge_chains(model, agchains,"Ag")
native = merge_chains(native, agchains, "Ag")

model = merge_chains(model, abchains, "Ab")
native = merge_chains(native, abchains, "Ab")


# native:model chain map dictionary for two interfaces
chain_map = {"Ag": "Ag", "Ab": "Ab"}

# returns a dictionary containing the results and the total DockQ score
dockq_results, total_dockq_score = DockQ.DockQ.run_on_all_native_interfaces(
    model, native, chain_map=chain_map
)


# Initialize an empty dictionary to store statistics for each metric
metrics = ["DockQ", "LRMSD", "iRMSD", "fnat", "clashes"]
output_data = {"sample_id": sample_id}
precision = 2

if len(dockq_results) > 1:
    sys.exit("ERROR: Expected a single interface.")


# Loop through each metric and calculate statistics
for metric in metrics:
    if "AgAb" in dockq_results.keys():
        output_data["abag_" + metric.lower()] = round(dockq_results["AgAb"][metric], precision)
        output_data["abag_receptor"] = "Ag"
    elif "AbAg" in dockq_results.keys():
        output_data["abag_" + metric.lower()] = round(dockq_results["AbAg"][metric], precision)
        output_data["abag_receptor"] = "Ab"
    else:
        sys.exit("ERROR: Unexpected interface.")


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

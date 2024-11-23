from snakemake.script import snakemake
import torch
import pandas as pd
import pickle


def load_data_from_pkl(pkl_path: str):
    """Load data from pickle file."""
    with open(pkl_path, "rb") as p:
        data = pickle.load(p)
    return data


precision = 2
output_data = {
    "sample_id": snakemake.wildcards.sample_id,
}


af_data = load_data_from_pkl(snakemake.input.af_data)
ptm = torch.from_numpy(af_data["ptm"])
iptm = torch.from_numpy(af_data["iptm"])
ranking_confidence = torch.tensor(
    af_data["ranking_confidence"]
)  # for some reason 'ranking_confidence' is not a numpy array (like ptm and iptm) but np.float type
num_recycles=torch.tensor(
    af_data["num_recycles"]
)


output_data.update(
    {
        "ptm": round(ptm.item(), precision),
        "iptm": round(iptm.item(), precision),
        "ranking_confidence": round(ranking_confidence.item(), precision),
        "num_recycles": num_recycles.item()
    }
)



# Save the output data to a CSV file
output_df = pd.DataFrame([output_data])
output_df.to_csv(snakemake.output[0], index=False)

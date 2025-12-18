import argparse
import sys
import pandas as pd

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Run DockQ analysis with merged chains.")
    parser.add_argument("--sample_id", required=True, help="Sample ID for processing.")
    parser.add_argument("--abag_ipsae", required=True, help="")
    parser.add_argument("--output_csv", required=True, help="Path to save the output CSV file.")
    return parser.parse_args()





def main():
    args = parse_arguments()

    # Read input sample data
    sample_id = args.sample_id


    df = pd.read_csv(args.abag_ipsae, delim_whitespace=True)

    # --- Select row where Type == "max" ---
    row_max = df[df["Type"] == "max"].iloc[0].to_dict()

    output_data = {"sample_id": sample_id}
    precision = 2

    # --- Extract max metrics ---
    for key in ["ipSAE", "ipSAE_d0chn", "ipSAE_d0dom", "ipTM_d0chn", "pDockQ", "pDockQ2", "LIS"]:
        if key in row_max:
            val = row_max[key]
            if isinstance(val, (int, float)):
                val = round(val, precision)
            output_data[f"abag_max_{key}"] = val

    # --- Extract asym metrics ---
    asym_rows = df[df["Type"] == "asym"].reset_index(drop=True)
    if len(asym_rows) >= 2:
        for i in range(2):
            asym = asym_rows.iloc[i].to_dict()
            for key in ["Chn1", "ipSAE", "ipSAE_d0chn", "ipSAE_d0dom", "ipTM_d0chn", "pDockQ", "pDockQ2", "LIS"]:
                if key in asym:
                    val = asym[key]
                    if isinstance(val, (int, float)):
                        val = round(val, precision)
                    output_data[f"abag_asym_{i+1}_{key}"] = val
    else:
        print("Warning: fewer than 2 asym rows found!")

    # Save the output data to a CSV file
    output_df = pd.DataFrame([output_data])
    output_df.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    main()

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

    # Select the row where Type == "max"
    row = df[df["Type"] == "max"]

    # Convert that single row to a dict
    row_dict = row.iloc[0].to_dict()

    output_data = {}
    output_data["sample_id"] = sample_id

    precision = 2

    for key, value in row_dict.items():
        if key in ["ipSAE", "ipSAE_d0chn" ,"ipSAE_d0dom",  "ipTM_d0chn"   ,  "pDockQ" ,    "pDockQ2"   , "LIS"]:
            output_data[f"abag_{key}"] = round(value, precision)


    # Save the output data to a CSV file
    output_df = pd.DataFrame([output_data])
    output_df.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    main()

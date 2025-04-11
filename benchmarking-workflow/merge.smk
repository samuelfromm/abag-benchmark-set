

import csv
import os
import pandas as pd
import functools
import glob


output_dir = config["output_dir"]
samples_csv = config["samples_csv"]


# Define the rule all to run all the outputs
rule all:
    input:
        os.path.join(output_dir,"scores_success.csv"),


# Merge all scores that exist
rule merge_successful_jobs:
    output:
        temp(output_dir+"/scores_prel_success.csv")
    run:
        csv_files = glob.glob(os.path.join(output_dir,"*","*_merged.csv"))

        dataframes = [pd.read_csv(f) for f in csv_files]
    

        merged_df = pd.concat(dataframes, ignore_index=True)
        merged_df = merged_df.sort_values(by='sample_id')
        
        # Save the output data to a CSV file
        merged_df.to_csv(output[0], index=False)


rule merge_final:
    input:
        samples_csv,
        output_dir+"/scores_prel_success.csv",
    output:
        output_dir+"/scores_success.csv"
    run:
        import pandas as pd
        import functools

        dataframes = [pd.read_csv(input[i]) for i in range(len(input))]

        # Merge dataframes on column 'sample_id'
        merged_df = functools.reduce(
            lambda left_df, right_df: pd.merge(left_df, right_df, on="sample_id", how="outer"),
            dataframes,
        )
        
        # Save the output data to a CSV file
        merged_df.to_csv(output[0], index=False)





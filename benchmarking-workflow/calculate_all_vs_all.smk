
'''

snakemake -s calculate_scores.smk --use-conda --configfile config/config_prebuilt.yaml 

NOTE: We run shell scripts instead of python scripts directly due to an issue with conda and snakemake
(see https://stackoverflow.com/questions/74479965/snakemake-doesnt-activate-conda-environment-correctly)
This does not happen on all systems that I tried but I was unable to locate the issue.
'''

import csv
import pandas as pd

output_dir = config["output_dir"]
samples_csv = config["samples_csv"]

input_samples = pd.read_csv(samples_csv, index_col="sample_id", na_values=[], keep_default_na=False)
samples = input_samples.index.tolist()


def get_sample_value(wildcards, column):
    return input_samples.loc[wildcards.sample_id, column]

# Define the rule all to run all the outputs
rule all:
    input:
        expand(output_dir+"/{sample_id}/{sample_id}_merged.csv", sample_id=input_samples.index),
        output_dir+"/scores.csv",


# Define the rule to run MM-align
rule run_get_chain_alignment:
    input:
        reference_pdb=lambda wildcards: get_sample_value(wildcards, 'reference_pdb'),
        query_pdb=lambda wildcards: get_sample_value(wildcards, 'query_pdb'),
    output:
        temp(output_dir+"/{sample_id}/{sample_id}_mmalign.csv")
    shell:
        """
        echo "CMD: {config[mmalign_path]} {input.reference_pdb} {input.query_pdb} -outfmt 2"
        outstr=$(timeout 5m {config[mmalign_path]} {input.reference_pdb} {input.query_pdb} -outfmt 2)

        # EXAMPLE OUTPUT MMALIGN
        # #PDBchain1	PDBchain2	TM1	TM2	RMSD	ID1	ID2	IDali	L1	L2	Lali
        # /home/sfromm/git/ae-dockq/data/tests/pdb_files/5nzz_1.pdb:A:B	/home/sfromm/git/ae-dockq/data/tests/pdb_files/5nzz_2.pdb:A:B	0.5164	0.5164	5.14	0.442	0.442	0.742	864	864	515
        # #Total CPU time is  2.30 seconds

        # Extract PDBchain1 and PDBchain2 values
        PDBchain1=$(echo "$outstr" | sed -n 2p | awk '{{print $1}}' | cut -d':' -f2-)
        PDBchain2=$(echo "$outstr" | sed -n 2p | awk '{{print $2}}' | cut -d':' -f2-)
        TM1=$(echo "$outstr" | sed -n 2p | awk '{{print $3}}' | cut -d':' -f2-)
        TM2=$(echo "$outstr" | sed -n 2p | awk '{{print $4}}' | cut -d':' -f2-)
        RMSD=$(echo "$outstr" | sed -n 2p | awk '{{print $5}}' | cut -d':' -f2-)

        echo "sample_id,TM_normalized_reference,TM_normalized_query,RMSD,aln_reference_cut,aln_query_cut\n{wildcards.sample_id},${{TM1}},${{TM2}},${{RMSD}},${{PDBchain1}},${{PDBchain2}}" > {output[0]}
        """

# Define the rule to run ANTIBODY-ANTIGEN DockQ between samples
rule run_pair_dockq:
    input:
        query_abag_path=lambda wildcards: get_sample_value(wildcards, 'query_pdb_abag'),
        reference_abag_path=lambda wildcards: get_sample_value(wildcards, 'reference_pdb_abag'),
    output:
        temp(output_dir+"/{sample_id}/{sample_id}_abag_pair_dockq.csv"),
    conda:
        config["run_dockq_env"]
    shell:
        """
        python scripts/run_pair_dockq.py \
            --sample_id {wildcards.sample_id} \
            --query_abag_path {input.query_abag_path} \
            --reference_abag_path {input.reference_abag_path} \
            --output_csv {output[0]}
        """



# Define rule to run mean(ae1 - ae2)
# Define rule to run mean(interface ae1 - interface ae2)
rule run_pair_mae:
    input:
        query_aligned_error=lambda wildcards: get_sample_value(wildcards, 'query_ae'),
        reference_aligned_error=lambda wildcards: get_sample_value(wildcards, 'reference_ae'),
    output:
        temp(output_dir+"/{sample_id}/{sample_id}_pair_mae.csv"),
    conda:
        config["run_calculate_aligned_error_env"]
    shell:
        """
        python scripts/run_pair_mae.py \
            --sample_id {wildcards.sample_id} \
            --query_aligned_error {input.query_aligned_error} \
            --reference_aligned_error {input.reference_aligned_error} \
            --output_csv {output[0]}
        """




# Define rule to run mean(pae1 - pae2)
# Define rule to run mean(interface pae1 - interface pae2)
rule run_pair_mpae:
    input:
        query_af_data=lambda wildcards: get_sample_value(wildcards, 'query_af_data'),
        query_af_features=lambda wildcards: get_sample_value(wildcards, 'query_af_features'),
        reference_af_data=lambda wildcards: get_sample_value(wildcards, 'reference_af_data'),
        reference_af_features=lambda wildcards: get_sample_value(wildcards, 'reference_af_features'),
    output:
        temp(output_dir+"/{sample_id}/{sample_id}_pair_mpae.csv"),
    conda:
        config["run_get_pae_prediction_env"]
    shell:
        """
        python scripts/run_pair_mpae.py \
            --sample_id {wildcards.sample_id} \
            --af_data_query {input.query_af_data} \
            --af_features_query {input.query_af_features} \
            --af_data_reference {input.reference_af_data} \
            --af_features_reference {input.reference_af_features} \
            --output_csv {output[0]}
        """




rule merge_scores:
    input:
        output_dir+"/{sample_id}/{sample_id}_mmalign.csv",
        output_dir+"/{sample_id}/{sample_id}_abag_pair_dockq.csv",
        output_dir+"/{sample_id}/{sample_id}_pair_mpae.csv",
        output_dir+"/{sample_id}/{sample_id}_pair_mae.csv",
    output:
        output_dir+"/{sample_id}/{sample_id}_merged.csv",
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


rule aggregate_sample_scores:
    input:
        expand(output_dir+"/{sample_id}/{sample_id}_merged.csv", sample_id=input_samples.index)
    output:
        temp(output_dir+"/scores_prel.csv")
    run:
        import pandas as pd

        dataframes = [pd.read_csv(input[i]) for i in range(len(input))]

        # Merge dataframes row wise
        merged_df = pd.concat(dataframes, ignore_index=True)
        
        merged_df = merged_df.sort_values(by='sample_id')
        
        # Save the output data to a CSV file
        merged_df.to_csv(output[0], index=False)

    
rule merge_final:
    input:
        output_dir+"/scores_prel.csv",
        samples_csv,
    output:
        output_dir+"/scores.csv"
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





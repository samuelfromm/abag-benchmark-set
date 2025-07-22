
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
        temp(output_dir+"/{sample_id}/{sample_id}_alignments.csv")
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

        echo "sample_id,aln_reference,aln_query\n{wildcards.sample_id},${{PDBchain1}},${{PDBchain2}}" > {output[0]}
        """



# Define the rule to run align and cut
rule run_align_and_cut:
    input:
        alignment=output_dir+"/{sample_id}/{sample_id}_alignments.csv",
        reference_pdb=lambda wildcards: get_sample_value(wildcards, 'reference_pdb'),
        query_pdb=lambda wildcards: get_sample_value(wildcards, 'query_pdb'),
    # params:
    #     reference_id=lambda wildcards: get_sample_value(wildcards, 'reference_id'),
    #     query_id=lambda wildcards: get_sample_value(wildcards, 'query_id'),
    output:
        temp(output_dir+"/{sample_id}/{sample_id}_align_and_cut.csv"),
        reference_cut=output_dir+"/{sample_id}/{sample_id}_cut_reference.pdb",
        query_cut=output_dir+"/{sample_id}/{sample_id}_cut_query.pdb",
    conda:
        config["run_align_and_cut_env"]
    shell:
        """
        python scripts/run_align_and_cut.py \
            --alignment {input.alignment} \
            --reference_pdb {input.reference_pdb} \
            --query_pdb {input.query_pdb} \
            --sample_id {wildcards.sample_id} \
            --reference_cut {output.reference_cut} \
            --query_cut {output.query_cut} \
            --output {output[0]}
        """



# Define the rule to run MM-align
rule run_mmalign:
    input:
        reference_cut=output_dir+"/{sample_id}/{sample_id}_cut_reference.pdb",
        query_cut=output_dir+"/{sample_id}/{sample_id}_cut_query.pdb",
    output:
        temp(output_dir+"/{sample_id}/{sample_id}_mmscore.csv")
    shell:
        """
        echo "CMD: {config[mmalign_path]} {input.reference_cut} {input.query_cut} -outfmt 2"
        outstr=$(timeout 5m {config[mmalign_path]} {input.reference_cut} {input.query_cut} -outfmt 2)

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


# Define the rule to run DockQ
rule run_dockq:
    input:
        alignment=output_dir+"/{sample_id}/{sample_id}_mmscore.csv",
        reference_cut=output_dir+"/{sample_id}/{sample_id}_cut_reference.pdb",
        query_cut=output_dir+"/{sample_id}/{sample_id}_cut_query.pdb",
    output:
        temp(output_dir+"/{sample_id}/{sample_id}_dockq.csv"),
    conda:
        config["run_dockq_env"]
    shell:
        """
        python scripts/run_dockq.py \
            --alignment {input.alignment} \
            --query_cut {input.query_cut} \
            --reference_cut {input.reference_cut} \
            --sample_id {wildcards.sample_id} \
            --output {output[0]}
        """


# Define the rule to run ANTIBODY-ANTIGEN DockQ
rule run_abag_dockq:
    input:
        samples_csv,
        reference_cut=output_dir+"/{sample_id}/{sample_id}_cut_reference.pdb",
        query_cut=output_dir+"/{sample_id}/{sample_id}_cut_query.pdb",
    output:
        temp(output_dir+"/{sample_id}/{sample_id}_abag_dockq.csv"),
    conda:
        config["run_dockq_env"]
    shell:
        """
        python scripts/run_abag_dockq.py \
            --input_csv {input[0]} \
            --query_pdb {input.query_cut} \
            --reference_pdb {input.reference_cut} \
            --sample_id {wildcards.sample_id} \
            --output {output[0]}
        """

# Define the rule to run DockQ
rule run_pdockq2:
    input:
        query_pdb=lambda wildcards: get_sample_value(wildcards, 'query_pdb'),
        af_data=lambda wildcards: get_sample_value(wildcards, 'query_af_data'),
    output:
        temp(output_dir+"/{sample_id}/{sample_id}_pdockq2.csv"),
    conda:
        config["run_pdockq2_env"]
    shell:
        """
        python scripts/run_pdockq2.py \
            --sample_id {wildcards.sample_id} \
            --query_pdb {input.query_pdb} \
            --af_data {input.af_data} \
            --output_csv {output[0]}
        """

rule run_get_af_prediction:
    input:
        af_data=lambda wildcards: get_sample_value(wildcards, 'query_af_data'),
    output:
        temp(output_dir+"/{sample_id}/{sample_id}_af_prediction.csv"),
    conda:
        config["run_get_af_prediction_env"]
    shell:
        """
        python scripts/run_get_af_prediction.py \
        --sample_id {wildcards.sample_id} \
        --af_data {input.af_data} \
        --output_csv {output[0]}
        """



# Calculate aligned error
rule run_calculate_aligned_error:
    input:
        reference_cut=output_dir+"/{sample_id}/{sample_id}_cut_reference.pdb",
        query_cut=output_dir+"/{sample_id}/{sample_id}_cut_query.pdb",
    output:
        output_dir+"/{sample_id}/{sample_id}_aligned_error.pth",
    conda:
        config["run_calculate_aligned_error_env"]
    shell:
        """
        python scripts/run_calculate_aligned_error.py \
        --sample_id {wildcards.sample_id} \
        --reference_pdb {input.reference_cut} \
        --query_pdb {input.query_cut} \
        --output_data {output[0]}
        """


# Calculate aeTM
rule run_calculate_aetm:
    input:
        aligned_error_data=output_dir+"/{sample_id}/{sample_id}_aligned_error.pth",
        reference_cut=output_dir+"/{sample_id}/{sample_id}_cut_reference.pdb",
        query_cut=output_dir+"/{sample_id}/{sample_id}_cut_query.pdb",
    output:
        temp(output_dir+"/{sample_id}/{sample_id}_aetm.csv"),
    conda:
        config["run_calculate_aetm_env"]
    shell:
        """
        python scripts/run_calculate_aetm.py \
        --sample_id {wildcards.sample_id} \
        --aligned_error_data {input.aligned_error_data} \
        --reference_pdb {input.reference_cut} \
        --query_pdb {input.query_cut} \
        --output_csv {output[0]}
        """

rule run_calculate_pae_prediction:
    input:
        af_data=lambda wildcards: get_sample_value(wildcards, 'query_af_data'),
        af_features=lambda wildcards: get_sample_value(wildcards, 'query_af_features'),
    output:
        temp(output_dir+"/{sample_id}/{sample_id}_pae_prediction.csv"),
    conda:
        config["run_get_pae_prediction_env"]
    shell:
        """
        python scripts/run_calculate_pae_prediction.py \
        --sample_id {wildcards.sample_id} \
        --af_data {input.af_data} \
        --af_features {input.af_features} \
        --output_csv {output[0]}
        """


rule merge_scores:
    input:
        output_dir+"/{sample_id}/{sample_id}_alignments.csv",
        output_dir+"/{sample_id}/{sample_id}_align_and_cut.csv",
        output_dir+"/{sample_id}/{sample_id}_dockq.csv",
        output_dir+"/{sample_id}/{sample_id}_abag_dockq.csv",
        output_dir+"/{sample_id}/{sample_id}_mmscore.csv",
        output_dir+"/{sample_id}/{sample_id}_pdockq2.csv",
        output_dir+"/{sample_id}/{sample_id}_af_prediction.csv",
        output_dir+"/{sample_id}/{sample_id}_aetm.csv",
        output_dir+"/{sample_id}/{sample_id}_pae_prediction.csv",
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





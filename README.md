# Introduction

This repository contains code used for the study ["Evaluating Deep Learning Based Structure Prediction Methods on Antibody-Antigen Complexes"](https://www.biorxiv.org/content/10.1101/2025.07.11.662141v1).

# Usage

To reproduce the results from the study or to apply the scoring pipeline to your own models, follow the instructions in the [Workflow](#workflow) section below.

As part of the paper, the benchmark dataset (including MSAs) is provided in a separate repository. [https://gitlab.com/ElofssonLab/abag-benchmark-dataset](https://gitlab.com/ElofssonLab/abag-benchmark-dataset) 

If you want to create your own version of the dataset used in our study — starting from [SAbDab](https://doi.org/10.1093/nar/gkab1050) — with different parameters, please refer to [https://github.com/samuelfromm/AADaM-fork](https://github.com/samuelfromm/AADaM-fork).

For access to the remaining data, including the model files, etc., please refer to the links in the publication.

## Workflow

1. Deposit the model files under `data/models`.
2. Depending on the output format of your files, you may need to rename them, e.g., run `scripts/rename_af3_files.py`.
3. Create input files for the benchmarking workflow; see for example `scripts-tetralith/run_create_samples_per_id_af3.sh`.
4. *(Optional)* Verify that the input files were correctly generated using `scripts-tetralith/check_scores_exist.sh`.
5. Run the benchmarking workflow:
    - With aligned error analysis: `scripts-tetralith/calculate_scores_ae_per_id.sh`
    - Without aligned error analysis: `scripts-tetralith/calculate_scores_per_id.sh`
6. Aggregate the results into a single file by running `scripts-tetralith/aggregate_scores_per_id.sh`.

The Conda environment used to run the workflow is defined in `snakemake_env.yml`.

## Sampling analysis

To run the sampling analysis from the study, run `scripts-tetralith/run_sampling_ae.sh` or `scripts-tetralith/run_sampling.sh`.
Additional scripts that might be helpful to analyse the results can be found at `analysis`.
The Conda environment used to run the analysis is defined in `analysis_env.yml`.

# Repository Structure

``````
abag-benchmark-set
├── README.md
├── analysis/                # Analysis scripts and notebooks
├── benchmarking-workflow/   # Snakemake workflow for generating benchmark metrics
├── data/                    # Raw and processed input data. For the complete data please refer to the information in the article
├── modules/                 # Submodules
├── results/                 # Output results from sampling analysis
├── scripts/                 # General-purpose scripts
├── scripts-berzelius/       # Scripts specific to the GPU cluster
├── scripts-tetralith/       # Scripts specific to the CPU cluster
├── src/                     # Source code for models or processing
│   ├── aligned_error        # Source code for aligned error metrics
│   ├── confidence           # Source code for confidence metrics
│   ├── bioutils             # Utility functions
``````
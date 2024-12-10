import pandas as pd
import random
import argparse
import numpy as np


def analyze_pdb_data_with_steps(data, max_sample_size, step_size, iterations, columns_to_analyze, reference_column, preset):
    """
    Analyzes PDB data by sampling rows and calculating averages for specified columns.

    Parameters:
    - data: DataFrame containing the PDB data.
    - max_sample_size: Maximum number of rows to sample.
    - step_size: Increment step for sample size.
    - iterations: Number of iterations to perform per sample size.
    - columns_to_analyze: List of column names to analyze.
    - reference_column: Column used for reference values.
    - preset: String to tag the results with a preset name.

    Returns:
    - A DataFrame containing the results of the analysis.
    """
    results = []
    pdbids = data['pdbid'].unique()

    for i, pdbid in enumerate(pdbids):
        print(f"Processing pdbid: {pdbid} ({i+1}/{len(pdbids)})")

        pdb_data = data[data['pdbid'] == pdbid]

        if len(pdb_data) < max_sample_size:
            print(f"Skipping pdbid {pdbid} (rows: {len(pdb_data)} < max_sample_size: {max_sample_size})")
            continue

        for sample_size in range(1, max_sample_size + 1, step_size):
            sample_results = {'pdbid': pdbid, 'sample_size': sample_size, 'preset': preset}

            column_results = {col: [] for col in columns_to_analyze}
            random_results = []

            for _ in range(iterations):
                # Sample rows with replacement
                sample = pdb_data.sample(n=sample_size, replace=True, random_state=random.randint(0, 10000))

                # Calculate max-based reference for each column
                for column in columns_to_analyze:
                    max_value = sample[column].max()
                    max_row = sample[sample[column] == max_value].iloc[0]
                    column_results[column].append(max_row[reference_column])


            # Calculate the averages
            for column in columns_to_analyze:
                sample_results[f'avg_{column}'] = np.mean(column_results[column])

            results.append(sample_results)

    return pd.DataFrame(results)


def parse_arguments():
    """
    Parses command-line arguments and returns them as a namespace object.
    """
    parser = argparse.ArgumentParser(description="Analyze PDB data with varying sample sizes.")
    parser.add_argument("--input_csv_paths", nargs="+", required=True, help="Paths to the input CSV files (space-separated)")
    parser.add_argument("--output_csv_path", required=True, help="Path to save the output CSV file")
    parser.add_argument("--max_sample_size", type=int, required=True, help="Maximum sample size")
    parser.add_argument("--step_size", type=int, default=1, help="Step size for sample size increments")
    parser.add_argument("--iterations", type=int, default=50, help="Number of iterations for averaging")
    parser.add_argument("--columns_to_analyze", nargs="+", required=True, help="Columns to analyze (space-separated)")
    parser.add_argument("--reference_column", required=True, help="Column for reference values (e.g., abag_dockq)")
    parser.add_argument("--preset", required=True, help="Preset name to include in the results")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    # Load and concatenate data from multiple CSV files
    data_frames = []
    for input_csv_path in args.input_csv_paths:
        try:
            df = pd.read_csv(input_csv_path)
            data_frames.append(df)
        except FileNotFoundError:
            print(f"File not found: {input_csv_path}")
    
    if not data_frames:
        raise ValueError("No valid CSV files were provided.")

    data = pd.concat(data_frames, ignore_index=True)

    # Analyze PDB data
    results_df = analyze_pdb_data_with_steps(
        data=data,
        max_sample_size=args.max_sample_size,
        step_size=args.step_size,
        iterations=args.iterations,
        columns_to_analyze=args.columns_to_analyze,
        reference_column=args.reference_column,
        preset=args.preset
    )

    # Save results
    results_df.to_csv(args.output_csv_path, index=False)
    print(f"Results saved to {args.output_csv_path}")
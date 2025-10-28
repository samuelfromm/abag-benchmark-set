import argparse
import sys
import os
import json
import pickle
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser, Chain, MMCIFIO


def load_structure(file_path, structure_id="structure"):
    """Load a PDB or mmCIF structure using Biopython."""
    if file_path.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(structure_id, file_path)


def merge_chains_dict(structure, merge_map):
    """
    Merge chains based on a dictionary mapping of new_chain_id -> list of input chain IDs.

    Example:
        merge_map = {
            "A": ["H", "L"],    # Merge H + L -> A
            "B": ["C", "D", "E"]  # Merge C + D + E -> B
        }
    """
    model = next(structure.get_models())  # assume single model
    new_chains = {}

    for new_chain_id, input_chain_ids in merge_map.items():
        new_chain = Chain.Chain(new_chain_id)
        residue_counter = 1

        for chain_id in input_chain_ids:
            try:
                chain = model[chain_id]
            except KeyError:
                print(f"Warning: Chain {chain_id} not found, skipping.")
                continue

            for res in chain:
                res.id = (' ', residue_counter, ' ')
                new_chain.add(res)
                residue_counter += 1

        new_chains[new_chain_id] = new_chain

    # Remove all original chains
    for chain_id in list(model.child_dict.keys()):
        model.detach_child(chain_id)

    # Add new merged chains
    for chain in new_chains.values():
        model.add(chain)

    return structure


def save_structure(structure, output_path):
    """Save structure as mmCIF (.cif) file."""
    if not output_path.endswith(".cif"):
        output_path += ".cif"

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(output_path)
    print(f"Saved merged structure as CIF: {output_path}")


def load_data(file_path):
    """Load data from a pickle (.pkl) or JSON (.json) file."""
    _, file_extension = os.path.splitext(file_path)

    if file_extension == ".pkl":
        with open(file_path, "rb") as f:
            data = pickle.load(f)
    elif file_extension == ".json":
        with open(file_path, "r") as f:
            data = json.load(f)
    else:
        sys.exit(f"Unsupported file format: {file_extension}. Please use a .pkl or .json file.")

    if not isinstance(data, dict):
        sys.exit("Error: Loaded data is not a dictionary.")

    print(f"Loaded data from {file_path} with keys: {list(data.keys())}")
    return data


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Merge chains in a structure using alignment and sample data.")
    parser.add_argument("--alignment", required=True, help="Path to the alignment CSV file (contains aln_reference and aln_query columns).")
    parser.add_argument("--sample_id", required=True, help="Sample ID for processing.")
    parser.add_argument("--input_csv", required=True, help="Path to input CSV file with sample chain info (Hchain, Lchain, Achain).")
    parser.add_argument("--query_pdb", required=True, help="Path to the query PDB or CIF file.")
    parser.add_argument("--query_data", required=True, help="Path to data file associated with the query structure (.pkl or .json).")
    parser.add_argument("--output_path", required=True, help="Path to save merged structure (will be written as .cif).")
    return parser.parse_args()


def main():
    args = parse_arguments()

    # Load alignment data
    alignment_df = pd.read_csv(args.alignment, index_col="sample_id")

    # Extract alignment for this sample
    sample_id = args.sample_id
    alignment = [
        alignment_df.loc[sample_id, "aln_reference"],
        alignment_df.loc[sample_id, "aln_query"],
    ]

    native_chain_ids = alignment[0].split(":")
    model_chain_ids = alignment[1].split(":")

    if "" in native_chain_ids or "" in model_chain_ids:
        print("WARNING: Not all chains are aligned. Cleaning alignment list.")
        native_chain_ids = [c for c in native_chain_ids if c]
        model_chain_ids = [c for c in model_chain_ids if c]

    if len(native_chain_ids) != len(model_chain_ids):
        sys.exit("ERROR: Alignment has different number of chains between reference and query.")

    # Load sample info
    input_samples_df = pd.read_csv(args.input_csv, index_col="sample_id", na_values=[], keep_default_na=False)

    try:
        Hchain = input_samples_df.at[sample_id, "Hchain"]
        Achain = input_samples_df.at[sample_id, "Achain"]
        Lchain = input_samples_df.at[sample_id, "Lchain"]
    except KeyError as e:
        sys.exit(f"Missing column in input CSV: {e}")

    agchains = [c.strip() for c in Achain.split("|")]
    abchains = [Hchain]
    if Lchain != "NA":
        abchains.append(Lchain)

    assert all(c in native_chain_ids for c in agchains), "Not all Ag chains found in native alignment."
    assert all(c in native_chain_ids for c in abchains), "Not all Ab chains found in native alignment."

    print(f"Ag chains: {agchains}")
    print(f"Ab chains: {abchains}")

    # Map native → model chain IDs via alignment
    agchains = [model_chain_ids[native_chain_ids.index(c)] for c in agchains]
    abchains = [model_chain_ids[native_chain_ids.index(c)] for c in abchains]

    # Define merge mapping dict
    merge_map = {
        "A": agchains,   # Antigen -> A
        "B": abchains    # Antibody -> B
    }

    # Load structure
    structure = load_structure(args.query_pdb, "model")

    # Merge chains
    merged_structure = merge_chains_dict(structure, merge_map)

    # Save merged structure
    save_structure(merged_structure, args.output_path)

    # Load query data
    data = load_data(args.query_data)

    # TODO: Validate atom_chain_ids vs merged_structure atom count
    assert len(data["atom_chain_ids"]) == len(list(merged_structure.get_atoms()))

    print("Done merging structure and loading data!")


if __name__ == "__main__":
    main()
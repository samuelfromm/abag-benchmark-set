from Bio import PDB
import sys
import argparse


def filter_pdb_by_chains(input_pdb, output_pdb, chains_to_keep):
    # Initialize the PDB parser and structure I/O
    parser = PDB.PDBParser(QUIET=True)
    io = PDB.PDBIO()

    # Parse the structure
    structure = parser.get_structure("structure", input_pdb)

    # Create a list to hold the selected chains
    selected_chains = []

    # Loop through models and chains, adding only the chains we want
    for model in structure:
        for chain in model:
            if chain.id in chains_to_keep:
                selected_chains.append(chain)

    # Define a new structure with selected chains
    class ChainSelect(PDB.Select):
        def accept_chain(self, chain):
            return chain in selected_chains

    # Write the new structure to the output PDB file
    io.set_structure(structure)
    io.save(output_pdb, ChainSelect())


def main():
    parser = argparse.ArgumentParser(description="Filter chains from a PDB file.")
    parser.add_argument("input_pdb", help="Input PDB file")
    parser.add_argument("output_pdb", help="Output PDB file")
    parser.add_argument("chains", nargs="+", help="List of chains to keep")

    args = parser.parse_args()

    filter_pdb_by_chains(args.input_pdb, args.output_pdb, args.chains)


if __name__ == "__main__":
    main()

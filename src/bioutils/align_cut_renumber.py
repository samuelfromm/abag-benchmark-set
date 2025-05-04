import sys
import argparse
import Bio.PDB
from Bio import Align
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio import SeqIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.Selection import unfold_entities
import Bio.SeqUtils
import pandas as pd
import io
from Bio.PDB import MMCIFParser, MMCIFIO
import os

### Global variables  - START ###

# SOURCE:   https://proteopedia.org/wiki/index.php/Amino_Acids
d3to1 = {
    "ALA": "A",
    "ASX": "B",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PYL": "O",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "SEC": "U",
    "VAL": "V",
    "TRP": "W",
    "UNK": "X",
    "TYR": "Y",
    "GLX": "Z",
}


### Global variables - END ###


def align_cut_renumber_sequences(model_chain, native_chain, rename_chains=None):
    """Takes two pdb chains and returns the part of the respective chains where the residues overlap (see example in the code below.)
        Args:
            model_chain: pdb chain structure
            native_chain: pdb chain structure
            rename_chains: Rule to rename chains ("native", "model", no renaming)

        Returns:
            new_model_chain: new pdb chain structrue
            new_native_chain: new pdb chain structure
            fraction_aligned_residues_native: percentage of aligned residues in the native structure
            fraction_aligned_residues_model: precentage of aligned residues in the model structure

    Example:
    model aligned aa sequence:      "NFS-SYDYG"
    native aligned aa sequence:     "--SSSYD-G"
    output_model:                   "SSYDG"
    numbering_model:                 12345
    output_native:                  "SSYDG"
    numbering_native:                12345
    """

    model_residues = unfold_entities(model_chain, "R")
    native_residues = unfold_entities(native_chain, "R")

    model_aa_sequence = "".join([d3to1[residue.resname] for residue in model_residues])
    native_aa_sequence = "".join(
        [d3to1[residue.resname] for residue in native_residues]
    )

    try:
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"

        # aligner.match_score = 5
        aligner.mismatch_score = (
            -100000
        )  # We do not allow any mismatch (to ensure that the aligned parts are identical)
        # aligner.open_gap_score = -10
        # aligner.extend_gap_score = -1
        # aligner.target_end_gap_score = 0.0
        # aligner.query_end_gap_score = 0.0

        # Perform the alignment
        alignments = aligner.align(seqA=model_aa_sequence, seqB=native_aa_sequence)
        optimal_alignment = alignments[0]  # Get the first (optimal) alignment


        aligned_model_aa_sequence = optimal_alignment[0,:]
        aligned_native_aa_sequence = optimal_alignment[1,:]

    except Exception as e:
        print(f"Got the following error: {e}")
        print(f"Trying to run Bio.pairwise2 instead")

        # pairwise2 is depreciated but the above code does not work for older versions of biopython
        alignments = Bio.pairwise2.align.globalms(
            model_aa_sequence,
            native_aa_sequence,
            match=1,
            mismatch=-100000,
            open=0,
            extend=0,
        )
        optimal_alignment = alignments[0]
        aligned_model_aa_sequence = optimal_alignment[0]
        aligned_native_aa_sequence = optimal_alignment[1]

    # renumber native chain residues and cut model to fit native chain structure
    if rename_chains == "native":
        # Use the name of the ***native*** chain for both the new model and native chain
        new_model_chain = Chain(native_chain.id)
        new_native_chain = Chain(native_chain.id)
    elif rename_chains == "model":
        # Use the name of the ***model*** chain for both the new model and native chain
        new_model_chain = Chain(model_chain.id)
        new_native_chain = Chain(model_chain.id)
    else:
        # keep the name of the chains as is
        new_model_chain = Chain(model_chain.id)
        new_native_chain = Chain(native_chain.id)

    model_residues_index = 0
    native_residues_index = 0
    new_joint_residue_id = 1  # start numbering residues from 1

    assert len(aligned_model_aa_sequence) == len(
        aligned_native_aa_sequence
    ), "Length of aligned model sequence does not equal length of aligned native sequence"

    for ix in range(
        len(optimal_alignment[1])
    ):  # Note len(alignment[1]) is equal to len(alignment[0])
        if (
            aligned_native_aa_sequence[ix] == "-"
            and aligned_model_aa_sequence[ix] != "-"
        ):  # skip current residue in the model
            model_residues_index += 1
            continue
        elif (
            aligned_model_aa_sequence[ix] == "-"
            and aligned_native_aa_sequence[ix] != "-"
        ):  # skip current residue in native structure
            native_residues_index += 1
            continue
        elif (
            aligned_model_aa_sequence[ix] != "-"
            and aligned_native_aa_sequence[ix] != "-"
        ):
            ## model rewrite
            model_residues[model_residues_index].detach_parent()
            model_residues[model_residues_index].id = (" ", new_joint_residue_id, " ")
            new_model_chain.add(model_residues[model_residues_index])

            ## native rewrite
            native_residues[native_residues_index].detach_parent()
            native_residues[native_residues_index].id = (" ", new_joint_residue_id, " ")
            new_native_chain.add(native_residues[native_residues_index])

            model_residues_index += 1
            native_residues_index += 1
            new_joint_residue_id += 1

    # check length of alignment

    num_residues_native_chain = len(list(native_chain.get_residues()))
    num_residues_model_chain = len(list(model_chain.get_residues()))

    num_residues_new_chain = len(list(new_model_chain.get_residues()))

    fraction_aligned_residues_native = float(num_residues_new_chain) / float(
        num_residues_native_chain
    )

    fraction_aligned_residues_model = float(num_residues_new_chain) / float(
        num_residues_model_chain
    )

    # if min(fraction_aligned_residues_native, fraction_aligned_residues_model) < 0.80:
    #     print(
    #         f"WARNING: Low number of aligned residues for native (model) chain {native_chain.id} ({model_chain.id})."
    #     )
    #     # sys.exit()

    return (
        new_model_chain,
        new_native_chain,
        fraction_aligned_residues_native,
        fraction_aligned_residues_model,
    )


def remove_hetatm(structure):
    """Remove hetero atoms in the structure.

    Args:
        structure: biopython structure object

    Returns:
        Returns a biopython structure object.
    """

    model_to_remove = []
    chain_to_remove = []
    residue_to_remove = []

    ## remove residues
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] != " ":  # check if hetero flag is not empty
                    residue_to_remove.append((model.id, chain.id, res.id))

    for residue in residue_to_remove:
        structure[residue[0]][residue[1]].detach_child(residue[2])

    ## remove chains
    for model in structure:
        for chain in model:
            if len(chain) == 0:
                chain_to_remove.append((model.id, chain.id))

    for chain in chain_to_remove:
        structure[chain[0]].detach_child(chain[1])

    ## remove models
    for model in structure:
        if len(model) == 0:
            model_to_remove.append(model.id)

    for model in model_to_remove:
        structure.detach_child(model)

    return structure


def add_arguments(parser):
    parser.add_argument(
        "-N", "--native", nargs=1, type=str, default="", help="path to model pdb/cif file"
    )
    parser.add_argument(
        "-M", "--model", nargs=1, type=str, default="", help="path to model pdb/cif file"
    )
    parser.add_argument(
        "-n",
        "--native_out",
        nargs=1,
        type=str,
        default="",
        help="path to native output file",
    )
    parser.add_argument(
        "-m",
        "--model_out",
        nargs=1,
        type=str,
        default="",
        help="path to model output file",
    )
    parser.add_argument(
        "-A",
        "--alignment",
        nargs=2,
        type=str,
        default=None,
        help='native chains and model chains superimposed onto each other. Chains are seperated by ":". Native chains come first. Example: "A:B:C" "L:H:A"',
    )
    parser.add_argument(
        "--rename_chains_rule",
        nargs=1,
        type=str,
        default=None,
        help='rename the cut native and model chains according to: native structure "native", model structure "model". If left out keep chain names as is.',
    )


def align_and_cut(
    path_to_native,
    path_to_model,
    path_native_out,
    path_model_out,
    alignment,
    rename_chains_rule=None,
    return_string=False,
):
    """Takes two pdb chains and returns the part of the respective chains where the residues overlap (see example in the code below.)
    Args:
        path_to_native:
        path_to_model:
        path_native_out:
        path_model_out:
        alignment: An alignment as output by MMalign. Example: A:B::C E:F:G::
        rename_chains_rule: Rule to rename chains ("native","native", no renaming)
        return_string: Whether or not to return the cut pdb files as a string.

    Returns:

        new_native_pdb_str,
        new_model_pdb_str,
        fraction_aligned_residues_native,
        fraction_aligned_residues_model,
        min_fraction_aligned_chain_residues_native,
        min_fraction_aligned_chain_residues_model,

    """
    def load_structure(path):
        ext = os.path.splitext(path)[1].lower()

        if ext == ".pdb":
            parser = PDBParser(QUIET=True)
        elif ext == ".cif":
            parser = MMCIFParser(QUIET=True)
        else:
            raise ValueError(f"Unsupported file extension: {ext}")

        return parser.get_structure("", path)

    iopdb = PDBIO()


    model_structure = load_structure(path_to_model)
    model_structure = remove_hetatm(model_structure)
    model_model = model_structure[0]

    native_structure = load_structure(path_to_native)
    native_structure = remove_hetatm(native_structure)
    native_model = native_structure[0]

    new_model_structure = Structure("")
    new_model_model = Model(0)
    new_native_structure = Structure("")
    new_native_model = Model(0)

    print(
        f"Running align and cut with alignment: native (model): {alignment[0]} ({alignment[1]})"
    )

    if alignment is None:
        model_chain_ids = [chain.get_id() for chain in model_model]
        native_chain_ids = [chain.get_id() for chain in native_model]
    else:
        native_chain_ids = alignment[0].split(":")
        model_chain_ids = alignment[1].split(":")
    if "" in native_chain_ids or "" in model_chain_ids:
        # example: A:B::C E:F:G:: -> use alignment: ABC EFG
        print(
            "WARNING: Not all chains are aligned. Replacing alignment with basic alignment order."
        )
        native_chain_ids = [chain_id for chain_id in native_chain_ids if chain_id != ""]
        model_chain_ids = [chain_id for chain_id in model_chain_ids if chain_id != ""]

    if len(native_chain_ids) != len(model_chain_ids):
        sys.exit("ERROR: Alignment has different number chains.")

    min_fraction_aligned_chain_residues_native = []
    min_fraction_aligned_chain_residues_model = []
    for i in range(len(native_chain_ids)):
        native_chain = native_model[native_chain_ids[i]]
        model_chain = model_model[model_chain_ids[i]]
        (
            new_model_chain,
            new_native_chain,
            fraction_aligned_residues_native,
            fraction_aligned_residues_model,
        ) = align_cut_renumber_sequences(model_chain, native_chain, rename_chains_rule)
        new_model_model.add(new_model_chain)
        new_native_model.add(new_native_chain)
        min_fraction_aligned_chain_residues_native.append(
            fraction_aligned_residues_native
        )
        min_fraction_aligned_chain_residues_model.append(
            fraction_aligned_residues_model
        )

    min_fraction_aligned_chain_residues_native = min(
        min_fraction_aligned_chain_residues_native
    )
    min_fraction_aligned_chain_residues_model = min(
        min_fraction_aligned_chain_residues_model
    )

    new_model_pdb_str = ""
    new_model_structure.add(new_model_model)
    iopdb.set_structure(new_model_structure)
    iopdb.save(path_model_out)
    if return_string:
        output = io.StringIO()
        iopdb.save(output)
        new_model_pdb_str = output.getvalue()
        output.close()

    new_native_pdb_str = ""
    new_native_structure.add(new_native_model)
    iopdb.set_structure(new_native_structure)
    iopdb.save(path_native_out)
    if return_string:
        output = io.StringIO()
        iopdb.save(output)
        new_native_pdb_str = output.getvalue()
        output.close()

    num_residues_native_structure = len(list(native_structure.get_residues()))
    num_residues_model_structure = len(list(model_structure.get_residues()))
    num_residues_new = len(list(new_native_structure.get_residues()))

    fraction_aligned_residues_native = float(num_residues_new) / float(
        num_residues_native_structure
    )

    fraction_aligned_residues_model = float(num_residues_new) / float(
        num_residues_model_structure
    )

    return (
        new_native_pdb_str,
        new_model_pdb_str,
        fraction_aligned_residues_native,
        fraction_aligned_residues_model,
        min_fraction_aligned_chain_residues_native,
        min_fraction_aligned_chain_residues_model,
    )


def main():

    parser = argparse.ArgumentParser(
        description="Takes two pdb files and cuts both structures (model 0) to match the intersection of their alignment. NOTE: If no chains are given the program will assume that the model and native chains are aligned. It will rename the model chains to match the native chains."
    )
    add_arguments(parser)
    args = vars(parser.parse_args(sys.argv[1:]))

    (
        _,
        _,
        fraction_aligned_residues_native,
        fraction_aligned_residues_model,
        min_fraction_aligned_chain_residues_native,
        min_fraction_aligned_chain_residues_model,
    ) = align_and_cut(
        path_to_native=args["native"][0],
        path_to_model=args["model"][0],
        path_native_out=args["native_out"][0],
        path_model_out=args["model_out"][0],
        alignment=args["alignment"],
        rename_chains_rule=args["rename_chains_rule"][0],
        return_string=False,
    )

    print(f"Native: {args['native'][0].split('/')[-1].split('.')[0]}")
    print(f"Model: {args['native'][0].split('/')[-1].split('.')[0]}")
    print(
        f"Fraction aligned residues in native structure: fraction_aligned_residues_native: {fraction_aligned_residues_native}"
    )
    print(
        f"Fraction aligned residues in model structure: fraction_aligned_residues_model: {fraction_aligned_residues_model}"
    )
    print(
        f"Minimum fraction aligned residues in native chains: min_perc_aligned_residues_native: {min_fraction_aligned_chain_residues_native}"
    )
    print(
        f"Minimum fraction aligned residues in model chains: min_perc_aligned_residues_model: {min_fraction_aligned_chain_residues_model}"
    )


if __name__ == "__main__":
    main()

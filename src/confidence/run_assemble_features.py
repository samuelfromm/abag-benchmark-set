import confidence_tools
from typing import Optional
import pickle
import os
import sys
import argparse
import numpy as np


def assemble_features(
    results_pkl_path: str,
    pdb_file_path: str,
    ligand: list,
    receptor: Optional[list] = None,
    output_results_pkl_path: Optional[str] = None,
    model_nber: Optional[int] = 0,
) -> dict:
    """Assemble new features for ligand and receptor pair.


    Args:



    Returns:

    """

    assert os.path.isfile(results_pkl_path) and os.path.isfile(pdb_file_path)
    results = confidence_tools.load_data_from_pkl(results_pkl_path)
    pdb_structure = confidence_tools.load_pdb_structure(pdb_file_path)

    asym_id = confidence_tools.compute_asym_id_from_pdb(pdb_structure, model_nber)

    # 1D features of shape [num_res]
    features_one_dim = {"plddt": results["plddt"]}

    # 2D features of shape [num_res, num_res] or features of shape [num_res, num_res, ???]
    features_two_dim = {
        "distances": confidence_tools.compute_distances_from_pdb(
            pdb_structure, model_nber
        ),
        "aligned_confidence_probs": results["aligned_confidence_probs"],
        "sequence_neighbour_adj": confidence_tools.get_sequence_neighbour(asym_id),
        "predicted_aligned_error": results["predicted_aligned_error"],
    }

    new_features = confidence_tools.create_ligand_receptor_features(
        asym_id=asym_id,
        ligand=ligand,
        receptor=receptor,
        features_one_dim=features_one_dim,
        features_two_dim=features_two_dim,
    )

    if not output_results_pkl_path is None:
        with open(output_results_pkl_path, "wb") as f:
            pickle.dump(new_features, f, protocol=4)

    return new_features


def add_arguments(parser):
    parser.add_argument(
        "--afmodel_pdb",
        help="path to af PDB model file",
        type=str,
    )
    parser.add_argument(
        "--results_pkl",
        help="path to results pkl file",
    )
    parser.add_argument(
        "--output_pkl",
        help="path, where to store the new features",
    )
    parser.add_argument(
        "--ligand",
        nargs="*",
        type=float,
        help="a list of chain indexes which make up the ligand",
    )
    parser.add_argument(
        "--receptor",
        nargs="*",
        type=float,
        help="a list of chain indexes which make up the receptor",
    )


def main():
    parser = argparse.ArgumentParser(
        description="Assemble new features for ligand receptor pair."
    )
    add_arguments(parser)
    args = parser.parse_args()

    assemble_features(
        results_pkl_path=args.results_pkl,
        pdb_file_path=args.afmodel_pdb,
        ligand=args.ligand,
        receptor=args.receptor,
        output_results_pkl_path=args.output_pkl,
    )


if __name__ == "__main__":
    main()

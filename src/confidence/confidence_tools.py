import Bio.PDB
import numpy as np
from typing import Optional, Tuple
import pickle
import os
import sys
import torch

def load_data_from_pkl(pkl_path: str):
    """Load data from pickle file."""
    with open(pkl_path, "rb") as p:
        data = pickle.load(p)
    return data


def load_pdb_structure(pdb_path: str) -> Bio.PDB.Structure.Structure:
    """Load pdb structure from path

    Args:
        pdb_path: Path to pdb file.

    Returns:
        pdb_structure: Returns the corresponding biopython structure object.
    """

    parser = Bio.PDB.PDBParser()
    pdb_structure = parser.get_structure("model", pdb_path)

    return pdb_structure


def compute_asym_id_from_pdb(
    pdb_structure: Bio.PDB.Structure.Structure, model_nber: Optional[int] = 0
) -> Tuple[np.ndarray, np.ndarray]:
    """Computes the asymmetric ID from the pdb model. Note that if the PDB file is missing residues, these residues will be missing in the asymmetric ID.

    Args:
        pdb_structure: Bio.PDB.Structure.Structure
        model_nber: model to be used from the pdb_structure (default 0)

    Returns:
        asym_id: [num_res], the assymetric ID
    """

    model = pdb_structure[model_nber]

    asym_id = []
    for i, chain in enumerate(model):
        for residue in chain:
            # print(residue.get_id()[1])
            asym_id.append(i + 1)

    return np.array(asym_id, dtype=np.float64)


def compute_distances_from_pdb(
    pdb_structure: Bio.PDB.Structure.Structure, model_nber: Optional[int] = 0
) -> Tuple[np.ndarray, np.ndarray]:
    # This uses the contact distance calculation from pDockQv2

    """Extracts the distances between residues from the pdb model.

    Args:
        pdb_structure: Bio.PDB.Structure.Structure
        model_nber: model to be used from the pdb_structure (default 0)

    Returns:
        distances: [num_res, num_res] a matrix containing distances between pairs of residues (in Å) (i.e. entry (i,j) contains the distance (in Å) between residue i and j)
    """

    model = pdb_structure[model_nber]

    c_alphas = []
    for i, chain in enumerate(model):
        for residue in chain:
            c_alphas.append(residue["CA"] if "CA" in residue else residue["CB"])

    num_res = len(c_alphas)
    distances = np.zeros(shape=(num_res, num_res))

    for i, ca_i in enumerate(c_alphas):
        for j, ca_j in enumerate(c_alphas):
            if j < i:
                distances[i, j] = abs(ca_i - ca_j)

    distances = distances + distances.T

    return distances


def compute_contact_mask(
    distances: np.ndarray,
    dist: Optional[float] = 8.0,
    self_contacts: Optional[bool] = False,
) -> np.ndarray:
    """Computes a mask containing all the pairs of residues that are in contact (i.e. the distance between them is below a certain threshold).

    Args:
      distances: [num_res, num_res] a matrix containing distances between pairs of residues (in Å) (i.e. entry (i,j) contains the distance (in Å) between residue i and j)
      dist: distance in Å used to define a contact (two residues are in contact if the distance between them is less or equal than dist)
      self_contacts: If False, a residue will be counted as not being in contact with itself. If True, use values from distances.


    Returns:
      contact_mask: [num_res, num_res] mask containing all pairs of residues with distances less or equal than dist
    """
    num_res = distances.shape[0]
    contact_mask = np.zeros(shape=(num_res, num_res), dtype=bool)
    contact_mask = distances <= dist

    if self_contacts == False:
        np.fill_diagonal(contact_mask, 0)

    return contact_mask


def compute_interface_mask(asym_id: np.ndarray) -> np.ndarray:
    """Computes a interface mask: entry (i,j) is 1 if i and j are residues from different chains, else entry (i,j) is zero.

    Args:
      asym_id: [num_res] the asymmetric unit ID - the chain ID.

    Returns:
      interface_mask: [num_res, num_res]
    """
    interface_mask = asym_id[:, None] != asym_id[None, :]

    return interface_mask


def compute_interface_contact_mask(
    contact_mask: np.ndarray, asym_id: np.ndarray
) -> np.ndarray:
    """Computes a mask containing all pairs of residues (i,j) that are in contact where residue i and j are from different chains.

    Args:
      contact_mask: [num_res, num_res] mask containing all pairs of residues that are in contact
      asym_id: [num_res] the asymmetric unit ID - the chain ID.

    Returns:
      interface_contact_mask: [num_res, num_res] the submask defined by the contacts in the interface region.
    """
    num_res = asym_id.shape[0]

    pair_mask = np.ones(shape=(num_res, num_res), dtype=bool)
    pair_mask *= asym_id[:, None] != asym_id[None, :]

    interface_contact_mask = np.multiply(contact_mask, pair_mask)

    return interface_contact_mask


def create_ligand_receptor_features(
    asym_id: np.array,
    ligand: list,
    receptor: Optional[list] = None,
    features_one_dim: Optional[dict] = {},
    features_two_dim: Optional[dict] = {},
) -> dict:
    """Takes a ligand and receptor and creates new features to match ligand and receptor. If no receptor is given, take the complement of the ligand (as a sorted list). For instances, if the asymmetric id is [1,1,1,2,2,3] and receptor is [1], ligand is [3] the new asymmetric unit will be [1,2,2,2] where 1 matches chain 3 and 2 matches chain 1. The features are cut up according to the new asymmetric unit ID. If ligand or receptor contains several chains, the new asymmetric unit ID will be assembled in that order.


    Args:
        asym_id: [num_res] the asymmetric unit ID - the chain ID.
        ligand: list , the index(es) of the chain(s) in the ligand.
        receptor: list , the index(es) of the chain(s) in the receptor (must be distinct from ligand).
        features_one_dim: A dictionary of 1D features of shape [num_res] that should be cut up.
        features_two_dim: A dictionary of 2D (or higher) features of shape [num_res, num_res (, ???)] that should be cut up.

    Returns:
        A dictionary containing the new asymmetric unit ID, the ligand and receptor list as well as the new features (infered from features_one_dim and features_two_dim).

    """

    num_res = asym_id.shape[0]

    chains = np.unique(asym_id).tolist()

    if receptor is None:
        receptor = [chain for chain in chains if not chain in set(ligand)]

    assert set(receptor) <= set(chains)
    assert set(ligand) <= set(chains)
    assert set(receptor).isdisjoint(set(ligand))

    # Assemble new asym_id containing two chains: ligand (1) and receptor (2)
    new_asym_id = []
    for idx in ligand:
        new_asym_id += np.count_nonzero(asym_id == idx) * [1]
    for idx in receptor:
        new_asym_id += np.count_nonzero(asym_id == idx) * [2]
    new_asym_id = np.array(new_asym_id, dtype=np.float64)
    new_num_res = new_asym_id.shape[0]

    # 1D features of shape [num_res]

    for key in features_one_dim.keys():
        assert features_one_dim[key].shape[0] == num_res
        new_feature = np.zeros(new_num_res)

        lower_idx = 0
        for idx in ligand + receptor:
            mask = np.isin(asym_id, idx)
            sub_feature = features_one_dim[key][mask]
            upper_idx = lower_idx + sub_feature.shape[0]
            new_feature[lower_idx:upper_idx] = sub_feature
            lower_idx = upper_idx
        assert new_feature.shape[0] == new_num_res
        features_one_dim[key] = new_feature

    # 2D features of shape [num_res, num_res] or features of shape [num_res, num_res, ???]

    for key in features_two_dim.keys():
        assert (
            features_two_dim[key].shape[0] == features_two_dim[key].shape[1] == num_res
        )

        new_shape = list(features_two_dim[key].shape)
        new_shape[0] = new_shape[1] = new_num_res
        new_shape = tuple(new_shape)
        new_feature = np.zeros(new_shape)

        lower_idx0 = 0
        for idx0 in ligand + receptor:
            lower_idx1 = 0
            upper_idx0 = (
                lower_idx0 + features_two_dim[key][np.isin(asym_id, idx0), :].shape[0]
            )
            for idx1 in ligand + receptor:
                sub_feature = features_two_dim[key][np.isin(asym_id, idx0), :][
                    :, np.isin(asym_id, idx1)
                ]
                upper_idx1 = lower_idx1 + sub_feature.shape[1]
                new_feature[lower_idx0:upper_idx0, lower_idx1:upper_idx1] = sub_feature
                lower_idx1 = upper_idx1
            lower_idx0 = upper_idx0
        assert new_feature.shape[0] == new_feature.shape[1] == new_num_res
        features_two_dim[key] = new_feature

    # Create feature dict
    new_features = {}
    new_features["ligand"] = ligand
    new_features["receptor"] = receptor
    new_features["asym_id"] = new_asym_id
    new_features.update(features_one_dim)
    new_features.update(features_two_dim)

    return new_features


def _calculate_bin_centers(breaks: np.ndarray):
    """Gets the bin centers from the bin edges.

    Args:
        breaks: [num_bins - 1] the error bin edges.

    Returns:
        bin_centers: [num_bins] the error bin centers.
    """
    step = breaks[1] - breaks[0]

    # Add half-step to get the center
    bin_centers = breaks + step / 2
    # Add a catch-all bin at the end.
    bin_centers = np.concatenate([bin_centers, [bin_centers[-1] + step]], axis=0)
    return bin_centers


def calculate_predicted_tm_score(
    aligned_confidence_probs: np.ndarray,
    asym_id: Optional[np.ndarray] = None,
    interface: bool = False,
    residue_weights: Optional[np.ndarray] = None,
    max_error_bin: Optional[float] = 31.0,
    num_bins: Optional[int] = 64,
    return_per_alignment=False,
    return_predicted_tm_term=False,
) -> np.ndarray:
    """Computes predicted TM alignment or predicted interface TM alignment. The function behaves the same way as the predicted_tm_score function in confidence.py but allows for the use of aligned_confidence_probs instead of logits.

    Args:
        aligned_confidence_probs: [num_res, num_res, num_bins] (see prediction_results['aligned_confidence_probs']) (computed internally using the output logits from
        PredictedAlignedErrorHead.)
        max_error_bin: part of the config of the AF model, see config file ['model']['heads']['predicted_aligned_error']['max_error_bin']
        num_bins: part of the config of the AF model, see config file ['model']['heads']['predicted_aligned_error']['num_bins']
        max_error_bin and num_bins are used to define the breaks
        residue_weights: [num_res] the per residue weights to use for the
        expectation.
        asym_id: [num_res] the asymmetric unit ID - the chain ID. Only needed for
        ipTM calculation, i.e. when interface=True.
        interface: If True, interface predicted TM score is computed.
        return_per_alignment: If True, return per_alignment terms instead of TM score
        return_predicted_tm_term: If True, return predicted_tm_term instead of TM score


    Returns:
        per_alignment: [num_res] The predicted TM alignments or the predicted iTM alignments.
        ptm_score: The predicted TM alignment or the predicted iTM score.
    """
    num_res = aligned_confidence_probs.shape[0]
    breaks = np.linspace(
        0.0, max_error_bin, num_bins - 1
    )  # see initialization of the PredictedAlignedErrorHead class in AlphaFold

    # residue_weights has to be in [0, 1], but can be floating-point, i.e. the
    # exp. resolved head's probability.
    if residue_weights is None:
        residue_weights = np.ones(num_res)

    bin_centers = _calculate_bin_centers(breaks)

    # Clip num_res to avoid negative/undefined d0.
    clipped_num_res = max(num_res, 19)

    # Compute d_0(num_res) as defined by TM-score, eqn. (5) in Yang & Skolnick
    # "Scoring function for automated assessment of protein structure template
    # quality", 2004: http://zhanglab.ccmb.med.umich.edu/papers/2004_3.pdf
    d0 = 1.24 * (clipped_num_res - 15) ** (1.0 / 3) - 1.8

    # Use passed probs
    probs = aligned_confidence_probs

    # TM-Score term for every bin.
    tm_per_bin = 1.0 / (1 + np.square(bin_centers) / np.square(d0))
    # E_distances tm(distance).
    predicted_tm_term = np.sum(probs * tm_per_bin, axis=-1)

    if return_predicted_tm_term:
        return predicted_tm_term

    pair_mask = np.ones(shape=(num_res, num_res), dtype=bool)
    if interface:
        pair_mask *= asym_id[:, None] != asym_id[None, :]

    predicted_tm_term *= pair_mask

    pair_residue_weights = pair_mask * (
        residue_weights[None, :] * residue_weights[:, None]
    )
    normed_residue_mask = pair_residue_weights / (
        1e-8 + np.sum(pair_residue_weights, axis=-1, keepdims=True)
    )
    per_alignment = np.sum(predicted_tm_term * normed_residue_mask, axis=-1)

    if return_per_alignment:
        return per_alignment

    return np.asarray(per_alignment[(per_alignment * residue_weights).argmax()])


def get_sequence_neighbour(asym_id: np.array):
    """Create sequence neighbour feature.

    Args:
        asym_id: [num_res] the asymmetric unit id

    Returns:
        neighbour_adj: [num_res,num_res,1] sequence neighbour feature (1 if residue i and j are adjecent in the amino acid sequence, else 0)
    """
    num_res = asym_id.shape[0]
    # Create sequence neighbour feature (for all chains)
    # edge attr (i,j) is 1 if residue i and j are "neighbours" in the amino acid chain and 0 otherwise
    if num_res == 1:
        neighbour_adj = np.array([[0]])
    else:
        neighbour_adj = np.diag(np.ones(num_res - 1), k=1) + np.diag(
            np.ones(num_res - 1), k=-1
        )
        for i in range(1, num_res):
            if asym_id[i] != asym_id[i - 1]:
                neighbour_adj[i, i - 1] = 0
                neighbour_adj[i - 1, i] = 0
    neighbour_adj = neighbour_adj[:, :, np.newaxis]

    return neighbour_adj


def calculate_predicted_aligned_error(
    aligned_confidence_probs: np.ndarray,
    max_error_bin: Optional[float] = 31.0,
    num_bins: Optional[int] = 64,
) -> Tuple[np.ndarray, np.ndarray]:
    """Calculates the predicted aligned error. Note: This function is usually not needed as the pae is part of the features.

    Args:
    aligned_confidence_probs: [num_res, num_res, num_bins] (see prediction_results['aligned_confidence_probs']) (computed internally using the output logits from
    PredictedAlignedErrorHead.)
    max_error_bin: part of the config of the AF model, see config file ['model']['heads']['predicted_aligned_error']['max_error_bin']
    num_bins: part of the config of the AF model, see config file ['model']['heads']['predicted_aligned_error']['num_bins']
    max_error_bin and num_bins are used to define the breaks

    Returns:
    predicted_aligned_error: [num_res, num_res] the expected aligned distance
        error for each pair of residues.
    """

    breaks = np.linspace(
        0.0, max_error_bin, num_bins - 1
    )  # see initialization of the PredictedAlignedErrorHead class in AlphaFold

    bin_centers = _calculate_bin_centers(breaks)

    # Tuple of expected aligned distance error and max possible error.
    return np.sum(aligned_confidence_probs * bin_centers, axis=-1)

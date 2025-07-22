import torch
import openfold.data.data_transforms
import openfold.np.protein
from openfold.utils.rigid_utils import Rigid
from typing import Optional
import numpy as np
import warnings

def calculate_asym_id(protein: openfold.np.protein.Protein) -> torch.Tensor:
    """Calculates the asymmetric unit id.
    Args:
        protein:

    Returns:
        asym_id: [num_res]  the asymmetric unit ID - the chain ID.
    """

    asym_id = None
    # Check whether monomer or multimer
    if np.unique(protein.chain_index).size > 1:
        asym_id = protein.chain_index + 1

    return torch.from_numpy(asym_id)


def calculate_aligned_error(
    reference_protein: openfold.np.protein.Protein,
    query_protein: openfold.np.protein.Protein,
    return_asym_id: False,
) -> torch.Tensor:
    """Calculates the aligned error between a query structure and a reference structure. The output is a non-symmetric matrix e_ij that captures the error in the position of the Calpha atom of residue j when the query and reference structures are aligned using the backbone frame of residue i.
    For more information on how exactly the aligned error and the backbone frames are computated, we refer to Section 1.8 in the supplementary material of [1].

    [1] Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. Nature 596, 583–589 (2021). https://doi.org/10.1038/s41586-021-03819-2

    Args:
        reference_protein:
        query_protein:

    Returns:
        aligned_error: [num_res, num_res] the aligned error between the query and reference structure.
        asym_id: [num_res]  the asymmetric unit ID - the chain ID.
    """

    def calculate_frames(protein_obj):
        protein = {
            key: torch.from_numpy(getattr(protein_obj, key.replace("all_", "")))
            for key in ["aatype", "all_atom_positions", "all_atom_mask"]
        }
        protein = openfold.data.data_transforms.atom37_to_frames(protein)
        return openfold.data.data_transforms.get_backbone_frames(protein)

    # Assertions for matching aatype and atom mask

    if not reference_protein.aatype.shape == query_protein.aatype.shape:
        raise ValueError(
            f"The amino-acid of the reference structure does not coincide with the amino-acid of the query structure. "
            f"Note that this might not be an issue depending on the use case, in which one may disable this error."
        )

    # aatype: np.ndarray  # [num_res]
    # Amino-acid type for each residue represented as an integer between 0 and
    # 20, where 20 is 'X'.
    if not (reference_protein.aatype == query_protein.aatype).all():
        raise ValueError(
            f"The amino-acid of the reference structure does not coincide with the amino-acid of the query structure. "
            f"Note that this might not be an issue depending on the use case, in which one may disable this error."
        )

    # atom_mask: np.ndarray  # [num_res, num_atom_type]
    # Binary float mask to indicate presence of a particular atom. 1.0 if an atom
    # is present and 0.0 if not. This should be used for loss masking.
    if not (reference_protein.atom_mask == query_protein.atom_mask).all():

        warnings.warn(
            f"Not all atoms present in the query structure are present in the reference structure or vise versa. "
        )

    if not (reference_protein.chain_index == query_protein.chain_index).all():
        raise ValueError(f"The chain indices do not coincide.")
    asym_id = None
    # Check whether monomer or multimer
    if np.unique(reference_protein.chain_index).size > 1:
        asym_id = reference_protein.chain_index + 1

    # Process both reference and query proteins
    reference_protein = calculate_frames(reference_protein)
    query_protein = calculate_frames(query_protein)

    # Compute aligned error
    aligned_error = calculate_aligned_error_from_tesor(
        query_backbone_rigid_tensor=query_protein["backbone_rigid_tensor"],
        reference_backbone_rigid_tensor=reference_protein["backbone_rigid_tensor"],
    )

    return aligned_error, torch.from_numpy(asym_id)


def calculate_aligned_error_from_tesor(
    query_backbone_rigid_tensor: torch.Tensor,
    reference_backbone_rigid_tensor: torch.Tensor,
) -> torch.Tensor:
    """Calculates the aligned error between two tensors.

    Args:
        query_backbone_rigid_tensor: Backbone tensor of the query structure.
        reference_backbone_rigid_tensor: Backbone tensor of the reference structure.

    Returns:
        [num_res, num_res] the aligned error between the query and reference structure.
    """

    query_affine = Rigid.from_tensor_4x4(query_backbone_rigid_tensor)
    reference_affine = Rigid.from_tensor_4x4(reference_backbone_rigid_tensor)

    def _points(affine):
        pts = affine.get_trans()[..., None, :, :]
        return affine.invert()[..., None].apply(pts)

    sq_diff = torch.sum(
        (_points(query_affine) - _points(reference_affine)) ** 2, dim=-1
    ).detach()

    return torch.sqrt(sq_diff)


def get_backbone_coordinates(protein: openfold.np.protein.Protein) -> torch.Tensor:
    """
    Extracts the backbone atom coordinates (Calpha and Cbeta) from a PDB string.

    If Calpha coordinates are unavailable for a residue, Cbeta coordinates are used.

    Returns:
        backbone_coordinates: A [num_residues, 3] numpy array of the backbone coordinates in Angstroms.
    """

    # Find the indices for Cα and Cβ atom types
    CA_index = openfold.np.residue_constants.atom_types.index("CA")
    CB_index = openfold.np.residue_constants.atom_types.index("CB")

    # Extract Cα and Cβ coordinates
    CA_coords = protein.atom_positions[:, CA_index, :]
    CB_coords = protein.atom_positions[:, CB_index, :]

    # Select Cα coordinates where available, otherwise select Cβ coordinates
    backbone_coordinates = np.where(
        np.all(CA_coords != 0, axis=-1, keepdims=True), CA_coords, CB_coords
    )

    return torch.from_numpy(backbone_coordinates)


# TODO update this with lddt function / plddt loss function calculations
def compute_distances(protein: openfold.np.protein.Protein) -> torch.Tensor:
    """
    Computes a distance matrix between residues based on their Calpha atoms,
    falling back to Cbeta if Calpha is not available.

    Returns:
        distances: A [num_residues, num_residues] matrix where each entry (i, j)
                   contains the distance (in Ångströms) between the Calpha (or Cbeta) atoms of residues i and j.
    """

    # Get backbone coordinates for residues
    backbone_coordinates = get_backbone_coordinates(protein)

    # Compute pairwise Euclidean distances between residues
    distances = torch.linalg.vector_norm(
        backbone_coordinates[:, None] - backbone_coordinates, dim=-1
    )

    return distances


def compute_contact_mask(
    distances: torch.Tensor,
    dist: Optional[float] = 8.0,
    self_contacts: Optional[bool] = False,
) -> torch.Tensor:
    """Computes a mask containing all the pairs of residues that are in contact (i.e. the distance between them is below a certain threshold).

    Args:
      distances: [num_res, num_res] a matrix containing distances between pairs of residues (in Å) (i.e. entry (i,j) contains the distance (in Å) between residue i and j)
      dist: distance in Å used to define a contact (two residues are in contact if the distance between them is less or equal than dist)
      self_contacts: If False, a residue will be counted as not being in contact with itself. If True, use values from distances.


    Returns:
      contact_mask: [num_res, num_res] mask containing all pairs of residues with distances less or equal than dist
    """
    contact_mask = distances <= dist

    if not self_contacts:
        contact_mask.fill_diagonal_(0)
    return contact_mask


def compute_interface_contact_mask(
    contact_mask: torch.Tensor, asym_id: torch.Tensor
) -> torch.Tensor:
    """Computes a mask containing all pairs of residues (i,j) that are in contact where residue i and j are from different chains.

    Args:
      contact_mask: [num_res, num_res] mask containing all pairs of residues that are in contact
      asym_id: [num_res] the asymmetric unit ID - the chain ID.

    Returns:
      interface_contact_mask: [num_res, num_res] the submask defined by the contacts in the interface region.
    """

    num_res = asym_id.shape[0]

    pair_mask = torch.ones((num_res, num_res), dtype=torch.bool)
    pair_mask *= asym_id[:, None] != asym_id[None, :]

    interface_contact_mask = contact_mask * pair_mask

    return interface_contact_mask


def compute_interface_mask(asym_id: torch.Tensor) -> torch.Tensor:
    """Computes a interface mask: entry (i,j) is 1 if i and j are residues from different chains, else entry (i,j) is zero.

    Args:
      asym_id: [num_res] the asymmetric unit ID - the chain ID.

    Returns:
      interface_mask: [num_res, num_res]
    """
    interface_mask = asym_id[:, None] != asym_id[None, :]

    return interface_mask


def compute_interface_contact_mask_from_protein(
    protein: openfold.np.protein.Protein,
) -> torch.Tensor:
    distances = compute_distances(
        protein=protein,
    )
    contact_mask = compute_contact_mask(distances=distances)
    asym_id = calculate_asym_id(protein)
    interface_contact_mask = compute_interface_contact_mask(
        contact_mask=contact_mask, asym_id=asym_id
    )

    return interface_contact_mask

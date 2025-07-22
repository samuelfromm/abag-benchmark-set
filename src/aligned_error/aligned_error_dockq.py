import torch
import openfold.np.protein
from openfold.utils.rigid_utils import Rigid
from typing import Optional


from openfold.np.residue_constants import atom_types


# NOTE
# FOR NOW WE ASSUME THE INPUT IS A DIMER FOR NOW


# There is no good reason for why I am using torch.Tensor instead of numpy arrays


def calculate_ae_dockq(
    aligned_error: torch.Tensor,
    asym_id: torch.Tensor,
    interface_contact_mask: torch.Tensor,
):
    """Computes a dockq score using the aligned error.

    Args:
      aligned_error: [num_res, num_res] the aligned error between the query and reference structure.
      interface_contact_mask: [num_res, num_res] the submask defined by the contacts in the interface region.
      asym_id: [num_res] the asymmetric unit ID - the chain ID.

    Returns:
      interface_contact_mask: [num_res, num_res] the submask defined by the contacts in the interface region.
    """

    ae_LRMS_scaled = calculate_ae_LRMS_scaled(
        aligned_error=aligned_error, asym_id=asym_id
    ).item()
    ae_iRMS_scaled = calculate_ae_iRMS_scaled(
        aligned_error=aligned_error, interface_contact_mask=interface_contact_mask
    ).item()

    ae_dockq = (ae_LRMS_scaled + ae_iRMS_scaled) / 2.0  # we ignore pFnat for now

    return ae_dockq, ae_LRMS_scaled, ae_iRMS_scaled


def calculate_ae_iRMS_scaled(
    aligned_error: torch.Tensor,
    interface_contact_mask: torch.Tensor,
    eps: float = 1e-8,
):
    if not aligned_error.shape == interface_contact_mask.shape:
        raise ValueError(
            f"aligned_error.shape ({aligned_error.shape}) != interface_contact_mask.shape ({interface_contact_mask.shape})"
        )

    num_res = aligned_error.shape[0]

    residue_weights = aligned_error.new_ones(num_res)

    assert torch.sum(residue_weights) == num_res

    d2 = 1.5

    predicted_tm_term = 1.0 / (1 + (aligned_error**2) / (d2**2))

    assert num_res == residue_weights.shape[-1]

    pair_mask = interface_contact_mask

    predicted_tm_term *= pair_mask

    pair_residue_weights = pair_mask * (
        residue_weights[None, :] * residue_weights[:, None]
    )
    denom = eps + torch.sum(pair_residue_weights, dim=-1, keepdims=True)
    normed_residue_mask = pair_residue_weights / denom
    per_alignment = torch.sum(predicted_tm_term * normed_residue_mask, dim=-1)

    weighted = per_alignment * residue_weights

    argmax = (weighted == torch.max(weighted)).nonzero()[0]
    return per_alignment[tuple(argmax)]


def calculate_ae_LRMS_scaled(
    aligned_error: torch.Tensor,
    asym_id: torch.Tensor,
    eps: float = 1e-8,
):
    if not aligned_error.shape[0] == asym_id.shape[0]:
        raise ValueError(
            f"aligned_error.shape[0] ({aligned_error.shape[0]}) != asym_id.shape[0] ({asym_id.shape[0]})"
        )

    interface = True

    num_res = aligned_error.shape[0]

    residue_weights = aligned_error.new_ones(num_res)

    assert torch.sum(residue_weights) == num_res

    d1 = 8.5

    predicted_tm_term = 1.0 / (1 + (aligned_error**2) / (d1**2))

    assert num_res == residue_weights.shape[-1]

    pair_mask = residue_weights.new_ones((num_res, num_res), dtype=torch.int32)
    if interface and (asym_id is not None):
        pair_mask *= (asym_id[:, None] != asym_id[None, :]).to(dtype=pair_mask.dtype)

    predicted_tm_term *= pair_mask

    pair_residue_weights = pair_mask * (
        residue_weights[None, :] * residue_weights[:, None]
    )
    denom = eps + torch.sum(pair_residue_weights, dim=-1, keepdims=True)
    normed_residue_mask = pair_residue_weights / denom
    per_alignment = torch.sum(predicted_tm_term * normed_residue_mask, dim=-1)

    per_chain_itm = {}
    for chain_idx in torch.unique(asym_id):
        chain_idx = chain_idx.item()

        chain_per_alignment = per_alignment.clone()
        chain_id_mask = asym_id == chain_idx
        chain_per_alignment[~chain_id_mask] = 0

        weighted = chain_per_alignment * residue_weights

        argmax = (weighted == torch.max(weighted)).nonzero()[0]
        per_chain_itm[chain_idx] = per_alignment[tuple(argmax)]

    unique_elements, counts = torch.unique(asym_id, return_counts=True)

    max_count_index = torch.argmax(counts)

    key_with_max_count = unique_elements[max_count_index].item()

    return per_chain_itm[key_with_max_count]

import torch
from typing import Optional


def compute_aetm(
    aligned_error: torch.Tensor,
    residue_weights: Optional[torch.Tensor] = None,
    asym_id: Optional[torch.Tensor] = None,
    interface: bool = False,
    eps: float = 1e-8,
) -> torch.Tensor:
    """Computes predicted TM alignment or predicted interface TM alignment score using the aligned error. Note that this score is not "predicted", but it is different than the TM score.

    Args:
    aligned_error: [num_res, num_res] the aligned error.
    residue_weights: [num_res] the per residue weights to use for the
        expectation.
    asym_id: [num_res] the asymmetric unit ID - the chain ID. Only needed for
        ipTM calculation, i.e. when interface=True.
    interface: If True, interface predicted TM score is computed.

    Returns:
    aetm_score: The predicted TM alignment or the predicted iTM score computed from the aligned error.
    """

    num_res = aligned_error.shape[0]

    if residue_weights is None:
        residue_weights = aligned_error.new_ones(num_res)

    clipped_num_res = max(torch.sum(residue_weights), 19)

    assert torch.sum(residue_weights) == num_res

    d0 = 1.24 * (clipped_num_res - 15) ** (1.0 / 3) - 1.8

    predicted_tm_term = 1.0 / (1 + (aligned_error**2) / (d0**2))

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

    weighted = per_alignment * residue_weights

    argmax = (weighted == torch.max(weighted)).nonzero()[0]
    return per_alignment[tuple(argmax)]


def calculate_ranking_confidence(ptm: torch.Tensor, iptm: torch.Tensor) -> torch.Tensor:
    assert ptm.dim() == iptm.dim() == 0
    return 0.2 * ptm + 0.8 * iptm


def compute_per_chain_aetm(
    aligned_error: torch.Tensor,
    asym_id: torch.Tensor,
    residue_weights: Optional[torch.Tensor] = None,
    eps: float = 1e-8,
) -> torch.Tensor:
    """

    Args:
    aligned_error: [num_res, num_res] the aligned error.
    residue_weights: [num_res] the per residue weights to use for the
        expectation.
    asym_id: [num_res] the asymmetric unit ID - the chain ID. Only needed for
        ipTM calculation, i.e. when interface=True.
    interface: If True, interface predicted TM score is computed.

    Returns:

    """

    num_res = aligned_error.shape[0]

    if residue_weights is None:
        residue_weights = aligned_error.new_ones(num_res)

    clipped_num_res = max(torch.sum(residue_weights), 19)

    assert torch.sum(residue_weights) == num_res

    d0 = 1.24 * (clipped_num_res - 15) ** (1.0 / 3) - 1.8

    predicted_tm_term = 1.0 / (1 + (aligned_error**2) / (d0**2))

    assert num_res == residue_weights.shape[-1]

    pair_mask = residue_weights.new_ones((num_res, num_res), dtype=torch.int32)

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
    return per_chain_itm



def compute_custom_aetm(
    aligned_error: torch.Tensor,
    residue_weights: Optional[torch.Tensor] = None,
    asym_id: Optional[torch.Tensor] = None,
    interface: bool = False,
    weights: Optional[torch.Tensor] = None,
    eps: float = 1e-8,
) -> torch.Tensor:
    """Computes predicted TM alignment or predicted interface TM alignment score using the aligned error. Note that this score is not "predicted", but it is different than the TM score.

    Args:
    aligned_error: [num_res, num_res] the aligned error.
    residue_weights: [num_res] the per residue weights to use for the
        expectation.
    asym_id: [num_res] the asymmetric unit ID - the chain ID. Only needed for
        ipTM calculation, i.e. when interface=True.
    interface: If True, interface predicted TM score is computed.
    weights:

    Returns:
    aetm_score: The predicted TM alignment or the predicted iTM score computed from the aligned error.
    """

    num_res = aligned_error.shape[0]

    if residue_weights is None:
        residue_weights = aligned_error.new_ones(num_res)

    clipped_num_res = max(torch.sum(residue_weights), 19)

    assert torch.sum(residue_weights) == num_res

    d0 = 1.24 * (clipped_num_res - 15) ** (1.0 / 3) - 1.8

    predicted_tm_term = 1.0 / (1 + (aligned_error**2) / (d0**2))

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
    per_alignment = per_alignment * weights


    weighted = per_alignment * residue_weights

    return weighted.mean()


def compute_per_chain_custom_aetm(
    aligned_error: torch.Tensor,
    asym_id: torch.Tensor,
    residue_weights: Optional[torch.Tensor] = None,
    weights: Optional[torch.Tensor] = None,
    eps: float = 1e-8,
) -> torch.Tensor:
    """

    Args:
    aligned_error: [num_res, num_res] the aligned error.
    residue_weights: [num_res] the per residue weights to use for the
        expectation.
    asym_id: [num_res] the asymmetric unit ID - the chain ID. Only needed for
        ipTM calculation, i.e. when interface=True.
    interface: If True, interface predicted TM score is computed.
    weights: 

    Returns:

    """

    num_res = aligned_error.shape[0]

    if residue_weights is None:
        residue_weights = aligned_error.new_ones(num_res)

    clipped_num_res = max(torch.sum(residue_weights), 19)

    assert torch.sum(residue_weights) == num_res

    d0 = 1.24 * (clipped_num_res - 15) ** (1.0 / 3) - 1.8

    predicted_tm_term = 1.0 / (1 + (aligned_error**2) / (d0**2))

    assert num_res == residue_weights.shape[-1]

    pair_mask = residue_weights.new_ones((num_res, num_res), dtype=torch.int32)

    pair_mask *= (asym_id[:, None] != asym_id[None, :]).to(dtype=pair_mask.dtype)

    predicted_tm_term *= pair_mask

    pair_residue_weights = pair_mask * (
        residue_weights[None, :] * residue_weights[:, None]
    )
    denom = eps + torch.sum(pair_residue_weights, dim=-1, keepdims=True)
    normed_residue_mask = pair_residue_weights / denom
    per_alignment = torch.sum(predicted_tm_term * normed_residue_mask, dim=-1)
    per_alignment = per_alignment * weights

    per_chain_itm = {}
    for chain_idx in torch.unique(asym_id):
        chain_idx = chain_idx.item()

        chain_per_alignment = per_alignment.clone()
        chain_id_mask = asym_id == chain_idx
        chain_per_alignment[~chain_id_mask] = 0

        weighted = chain_per_alignment * residue_weights

        mean_weighted = weighted.mean()

        per_chain_itm[chain_idx] = mean_weighted
    return per_chain_itm

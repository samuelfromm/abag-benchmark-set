import argparse
import torch
import pandas as pd
import aligned_error.aligned_error_ops
import aligned_error.aligned_error_utils
import aligned_error.aligned_error_confidence

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Run aligned error confidence analysis.")
    parser.add_argument("--sample_id", required=True, help="Sample ID for processing.")
    parser.add_argument("--query_aligned_error", required=True, help="Path to the query aligned error file.")
    parser.add_argument("--reference_aligned_error", required=True, help="Path to the reference aligned error file.")
    parser.add_argument("--output_csv", required=True, help="Path to save the output CSV file.")
    return parser.parse_args()




def reorder_ae(ae, asym_id, new_chain_order):
    """
    Reorders a matrix and remaps chain labels based on a desired chain order (PyTorch version).

    Parameters:
    - ae: (n x n) torch tensor
    - asym_id: list or 1D torch tensor of original chain labels
    - new_chain_order: desired chain label order (e.g., [3, 1, 2])

    Returns:
    - reordered_ae: matrix with permuted rows and columns
    - new_asym_id: new asym_id with relabeled chains
    """
    index_list = []
    new_asym_id = torch.zeros_like(asym_id)

    # Build new index order and relabel chains
    for new_label, original_chain in enumerate(new_chain_order, start=1):
        indices = torch.where(asym_id == original_chain)[0]
        index_list.append(indices)
        new_asym_id[indices] = new_label  # relabel to 1, 2, 3...

    # Concatenate index groups into a single index tensor
    new_order = torch.cat(index_list)

    # Permute rows and columns
    reordered_ae = ae[new_order][:, new_order]

    # Reorder chain labels to match new matrix order
    new_asym_id = new_asym_id[new_order]

    return reordered_ae, new_asym_id



def main():
    args = parse_arguments()
    
    query_aligned_error = torch.load(args.query_aligned_error)
    query_ae = query_aligned_error["aligned_error"]
    query_asym_id = query_aligned_error["asym_id"]


    reference_aligned_error = torch.load(args.reference_aligned_error)
    reference_ae = reference_aligned_error["aligned_error"]
    reference_asym_id = reference_aligned_error["asym_id"]

    if not query_ae.shape == reference_ae.shape or not query_asym_id.shape == reference_asym_id.shape:
        output_data = {
            "sample_id": args.sample_id,
            "mae_interface": float("nan"),
            "mae": float("nan"),
            }
            # Save the output data to a CSV file
        output_df = pd.DataFrame([output_data])
        output_df.to_csv(args.output_csv, index=False)
        return
    
    if not torch.equal(reference_asym_id, query_asym_id):
        new_chain_order = reference_asym_id.unique().tolist()
        query_ae, query_asym_id = reorder_ae(query_ae, query_asym_id, new_chain_order)
        reference_ae, reference_asym_id = reorder_ae(reference_ae, reference_asym_id, new_chain_order)


    interface_contact_mask = aligned_error.aligned_error_ops.compute_interface_mask(
        asym_id=query_asym_id
    )


    mae_interface = (torch.abs(query_ae - reference_ae))[interface_contact_mask == 1].mean()
    mae = torch.abs(query_ae - reference_ae).mean()



    precision = 2
    # Prepare the output data
    output_data = {
        "sample_id": args.sample_id,
        "mae_interface": round(mae_interface.item(), precision),
        "mae": round(mae.item(), precision),
    }

    # Save the output data to a CSV file
    output_df = pd.DataFrame([output_data])
    output_df.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    main()

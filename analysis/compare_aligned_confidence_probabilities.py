import argparse
import confidence_tools

# import confidence
# import custom_confidence
import numpy as np
from typing import Dict, Optional, Tuple
import os, sys, inspect

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits import mplot3d

from copy import copy

import string

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

# import MMalign_wrapper
# import dgram2dmap_mod

import logging

logging.basicConfig(format="%(message)s", level=logging.INFO)


def add_arguments(parser):
    parser.add_argument(
        "--pkl1",
        help="path to pkl file for first model",
        type=str,
    )
    parser.add_argument(
        "--pkl2",
        help="path to pkl file for second model",
        type=str,
    )
    parser.add_argument(
        "--name1",
        help="name for first model",
        required=False,
        type=str,
        default="model1",
    )
    parser.add_argument(
        "--name2",
        help="name for second model",
        required=False,
        type=str,
        default="model2",
    )
    parser.add_argument(
        "--outfolder",
        help="outfolder",
        default=str,
    )
    parser.add_argument(
        "--MMalign_exe",
        help="path to MMalign executable",
        required=False,
        default=None,
    )


def get_chain_limits(asym_id):
    chain_limits = {}

    chain_ids = asym_id.astype("int")
    for i in range(chain_ids[-1]):
        chain_starts = np.where(chain_ids == i + 1)[0][0] + 1

        chain_stops = np.where(chain_ids == i + 1)[0][-1] + 1

        chain_limits[string.ascii_uppercase[i]] = (chain_starts, chain_stops)
    return chain_limits


def get_bounding_boxes(limitA, limitB, color="r"):
    rect1 = patches.Rectangle(
        (limitA[0] - 1, limitB[0] - 1),
        limitA[1] - limitA[0],
        limitB[1] - limitB[0],
        linewidth=1,
        edgecolor=color,
        facecolor="none",
    )
    rect2 = patches.Rectangle(
        (limitB[0] - 1, limitA[0] - 1),
        limitB[1] - limitB[0],
        limitA[1] - limitA[0],
        linewidth=1,
        edgecolor=color,
        facecolor="none",
    )
    return rect1, rect2


def plot_vec(
    filepath, vec1, vec2, limitA=None, limitB=None, name1="vec1", name2="vec2"
):
    l = vec1.shape[0]
    assert l == vec2.shape[0]
    difference = np.absolute(vec1 - vec2)

    # Figure Size
    fig, ax = plt.subplots(3, 1, figsize=(1 * 7, 2 * 7))
    X = np.arange(l)

    ax[0].bar(
        X,
        vec1 - vec2,
        color="b",
        width=1,
    )
    ax[0].legend(labels=[f"{name1}-{name2}"])

    if limitA is not None:

        ax00 = ax[0].twinx()
        # Plot a line
        ax00.axvline(x=limitA[1], color="r")

        ax00.annotate(
            f"{limitA[0]}-{limitA[1]}",
            xy=(limitA[1], 0),
            xytext=(limitA[1], -0.1),
            arrowprops=dict(facecolor="black", shrink=0.05),
        )

    ax[1].bar(X + 0.00, vec1, color="b", width=1)
    ax[1].bar(X + 0.0, vec2, color="g", width=1)
    ax[1].legend(labels=[name1, name2])

    if limitA is not None:
        ax11 = ax[1].twinx()
        # Plot a line
        ax11.axvline(x=limitA[1], color="r")

        ax11.annotate(
            f"{limitA[0]}-{limitA[1]}",
            xy=(limitA[1], 0),
            xytext=(limitA[1], -0.1),
            arrowprops=dict(facecolor="black", shrink=0.05),
        )

    ax[2].bar(X + 0.00, vec2, color="g", width=1)
    ax[2].bar(X + 0.0, vec1, color="b", width=1)
    ax[2].legend(labels=[name2, name1])

    if limitA is not None:
        ax22 = ax[2].twinx()
        # Plot a line
        ax22.axvline(x=limitA[1], color="r")

        ax22.annotate(
            f"{limitA[0]}-{limitA[1]}",
            xy=(limitA[1], 0),
            xytext=(limitA[1], -0.1),
            arrowprops=dict(facecolor="black", shrink=0.05),
        )

    plt.savefig(filepath, dpi=600)
    plt.close()


def plot_matrix(
    mat1, mat2, filepath=None, limitA=None, limitB=None, name1="mat1", name2="mat2"
):

    color_map = "inferno"
    absdifference = np.absolute(mat1 - mat2)
    difference1 = mat1 - mat2
    difference2 = mat2 - mat1

    fig, ax = plt.subplots(1, 5, figsize=(5 * 7, 1 * 7))
    p1 = ax[0].imshow(mat1, cmap=color_map)
    plt.colorbar(p1, ax=ax[0], fraction=0.046, pad=0.05)
    p2 = ax[1].imshow(mat2, cmap=color_map)
    plt.colorbar(p2, ax=ax[1], fraction=0.046, pad=0.05)
    ad = ax[2].imshow(absdifference, cmap=color_map)
    plt.colorbar(ad, ax=ax[2], fraction=0.046, pad=0.05)
    d1 = ax[3].imshow(difference1, cmap=color_map)
    plt.colorbar(d1, ax=ax[3], fraction=0.046, pad=0.05)
    d2 = ax[4].imshow(difference2, cmap=color_map)
    plt.colorbar(d2, ax=ax[4], fraction=0.046, pad=0.05)
    ax[0].title.set_text(name1)
    ax[1].title.set_text(name2)
    ax[2].title.set_text("absolute difference")
    ax[3].title.set_text(f"{name1} - {name2}")
    ax[4].title.set_text(f"{name2} - {name1}")

    if limitA and limitB:
        rect1, rect2 = get_bounding_boxes(limitA, limitB)
        rect3 = copy(rect1)
        rect4 = copy(rect2)
        rect5 = copy(rect1)
        rect6 = copy(rect2)
        rect7 = copy(rect1)
        rect8 = copy(rect2)
        rect9 = copy(rect1)
        rect10 = copy(rect2)

        ax[0].add_patch(rect1)
        ax[0].add_patch(rect2)
        ax[1].add_patch(rect3)
        ax[1].add_patch(rect4)
        ax[2].add_patch(rect5)
        ax[2].add_patch(rect6)
        ax[3].add_patch(rect7)
        ax[3].add_patch(rect8)
        ax[4].add_patch(rect9)
        ax[4].add_patch(rect10)

    if filepath is not None:
        plt.savefig(filepath, dpi=600)
        plt.close()
    else:
        return fig


def plot_matrix_3D(
    filepath, mat1, mat2, limitA=None, limitB=None, name1="mat1", name2="mat2"
):

    absdifference = np.absolute(mat1 - mat2)
    difference1 = mat1 - mat2
    difference2 = mat2 - mat1

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.grid()
    x = range(difference1.shape[0])
    y = range(difference1.shape[1])

    x, y = np.meshgrid(x, y)

    z = difference1 * 10
    ax.plot_surface(x, y, z[x, y])

    plt.savefig(filepath, dpi=600)
    plt.close()


def write_results(filepath, resultdict):
    keys = []
    values = []
    for key, value in resultdict.items():
        keys.append(key)
        values.append(str(value))

    output = ",".join(keys) + "\n" + ",".join(values)

    with open(filepath, "w") as file:
        file.write(output)


def plot_predicted_tm_term(
    aligned_confidence_probs_1: np.ndarray,
    aligned_confidence_probs_2: np.ndarray,
    asym_id_1: Optional[np.ndarray] = None,
    asym_id_2: Optional[np.ndarray] = None,
    filepath=None,
):

    aligned_confidence_probs_lst = [
        aligned_confidence_probs_1,
        aligned_confidence_probs_2,
    ]
    asym_id_lst = [asym_id_1, asym_id_2]

    predicted_tm_term_lst = [
        confidence_tools.calculate_predicted_tm_score(
            aligned_confidence_probs, return_predicted_tm_term=True
        )
        for aligned_confidence_probs in aligned_confidence_probs_lst
    ]

    if not all(asym_id is None for asym_id in asym_id_lst):
        chain_limits_lst = [get_chain_limits(asym_id) for asym_id in asym_id_lst]

    # Check if everything has the correct dimensions
    # assert chain_limits == chain_limits_2
    # check aligned_confidence_prbs have the same dimension

    limitA = chain_limits_lst[0]["A"]
    limitB = chain_limits_lst[0]["B"]

    fig = plot_matrix(
        predicted_tm_term_lst[0],
        predicted_tm_term_lst[1],
        name1=f"predicted_tm_term_1",
        name2=f"predicted_tm_term_2",
        filepath=filepath,
    )
    if filepath is None:
        return fig
    return None


def plot_predicted_aligned_error(
    aligned_confidence_probs_1: np.ndarray,
    aligned_confidence_probs_2: np.ndarray,
    asym_id_1: Optional[np.ndarray] = None,
    asym_id_2: Optional[np.ndarray] = None,
    filepath=None,
):

    aligned_confidence_probs_lst = [
        aligned_confidence_probs_1,
        aligned_confidence_probs_2,
    ]
    asym_id_lst = [asym_id_1, asym_id_2]

    predicted_tm_term_lst = [
        confidence_tools.calculate_predicted_tm_score(
            aligned_confidence_probs, return_predicted_tm_term=True
        )
        for aligned_confidence_probs in aligned_confidence_probs_lst
    ]

    if not all(asym_id is None for asym_id in asym_id_lst):
        chain_limits_lst = [get_chain_limits(asym_id) for asym_id in asym_id_lst]

    # Check if everything has the correct dimensions
    # assert chain_limits == chain_limits_2
    # check aligned_confidence_prbs have the same dimension

    limitA = chain_limits_lst[0]["A"]
    limitB = chain_limits_lst[0]["B"]

    fig = plot_matrix(
        predicted_tm_term_lst[0],
        predicted_tm_term_lst[1],
        name1=f"predicted_tm_term_1",
        name2=f"predicted_tm_term_2",
        filepath=filepath,
    )
    if filepath is None:
        return fig
    return None


# def main():

#     parser = argparse.ArgumentParser(description="Compare two models.")
#     add_arguments(parser)
#     args = parser.parse_args()

#     pickle_file_path_1 = args.pkl1
#     pickle_file_path_2 = args.pkl2

#     prediction_result_1 = confidence_tools.load_results_from_pkl(pickle_file_path_1)
#     prediction_result_2 = confidence_tools.load_results_from_pkl(pickle_file_path_2)

#     logging.info(f"Pkl 1: {os.path.split(pickle_file_path_1)[1]}")
#     logging.info(f"Pkl 2: {os.path.split(pickle_file_path_2)[1]}")

#     metrics = {}

#     # Compute scores
#     logging.info("Calculation alphafold scores.")
#     for i, prediction_result in enumerate([prediction_result_1, prediction_result_2]):
#         ptm = confidence.predicted_tm_score(
#             logits=prediction_result["predicted_aligned_error_raw"]["logits"],
#             breaks=prediction_result["predicted_aligned_error_raw"]["breaks"],
#         )
#         iptm = confidence.predicted_tm_score(
#             logits=prediction_result["predicted_aligned_error_raw"]["logits"],
#             breaks=prediction_result["predicted_aligned_error_raw"]["breaks"],
#             asym_id=prediction_result["predicted_aligned_error_raw"]["asym_id"],
#             interface=True,
#         )
#         ranking_confidence = 0.8 * iptm + 0.2 * ptm
#         metrics[f"ptm_{i}"] = ptm
#         metrics[f"iptm_{i}"] = iptm
#         metrics[f"ranking_confidence_{i}"] = ranking_confidence

#     # plot predicted aligned error and difference
#     logging.info("Calculating predicted tm term.")
#     predicted_tm_term_1 = custom_confidence.compute_predicted_tm_term(
#         logits=prediction_result_1["predicted_aligned_error_raw"]["logits"],
#         breaks=prediction_result_1["predicted_aligned_error_raw"]["breaks"],
#     )
#     predicted_tm_term_2 = custom_confidence.compute_predicted_tm_term(
#         logits=prediction_result_2["predicted_aligned_error_raw"]["logits"],
#         breaks=prediction_result_2["predicted_aligned_error_raw"]["breaks"],
#     )

#     logging.info("Computing chain limits.")
#     chain_limits = get_chain_limits(
#         prediction_result_1["predicted_aligned_error_raw"]["asym_id"]
#     )
#     chain_limits_2 = get_chain_limits(
#         prediction_result_2["predicted_aligned_error_raw"]["asym_id"]
#     )
#     logging.info(f"chain_limits_1: {chain_limits}")
#     logging.info(f"chain_limits_2: {chain_limits_2}")

#     assert chain_limits == chain_limits_2

#     limitA = chain_limits["A"]
#     limitB = chain_limits["B"]

#     logging.info("Plotting predicted tm term.")
#     plotfilepath = os.path.join(args.outfolder, "predicted_tm_term.png")

#     plot_matrix(
#         plotfilepath,
#         predicted_tm_term_1,
#         predicted_tm_term_2,
#         limitA,
#         limitB,
#         name1=f"predicted_tm_term_{args.name1}",
#         name2=f"predicted_tm_term_{args.name2}",
#     )

#     #####
#     pdb_file_path_1 = os.path.join(
#         os.path.split(pickle_file_path_1)[0],
#         "_".join(
#             ["unrelaxed"] + os.path.split(pickle_file_path_1)[1].split("_")[1:]
#         ).split(".")[0]
#         + ".pdb",
#     )
#     pdb_file_path_2 = os.path.join(
#         os.path.split(pickle_file_path_2)[0],
#         "_".join(
#             ["unrelaxed"] + os.path.split(pickle_file_path_2)[1].split("_")[1:]
#         ).split(".")[0]
#         + ".pdb",
#     )

#     if os.path.exists(pdb_file_path_1) and os.path.exists(pdb_file_path_2):
#         logging.info("Plotting predicted tm term with interface contact mask.")
#         """
#         # interface
#         asym_id = prediction_result_1['predicted_aligned_error_raw']['asym_id']
#         asym_id_2 = prediction_result_2['predicted_aligned_error_raw']['asym_id']
#         assert np.array_equal(asym_id, asym_id_2)
#         num_res = asym_id.shape[0]
#         pair_mask = np.ones(shape=(num_res, num_res), dtype=bool)
#         pair_mask *= asym_id[:, None] != asym_id[None, :]
#         predicted_tm_term_i1 = predicted_tm_term_1*pair_mask
#         predicted_tm_term_i2 = predicted_tm_term_2*pair_mask
#         plot_matrix(os.path.join(args.outfolder, 'predicted_tm_term_interface.png'),predicted_tm_term_i1,predicted_tm_term_i2,limitA,limitB,name1=f"predicted_tm_term_{args.name1}",name2=f"predicted_tm_term_{args.name2}")
#         """

#         pdb_structure_1 = confidence_tools.load_structure(pdb_file_path_1)
#         pdb_distances_1, pdb_asym_id_1 = confidence_tools.compute_distances_from_pdb(
#             pdb_structure_1
#         )

#         pdb_structure_2 = confidence_tools.load_structure(pdb_file_path_2)
#         pdb_distances_2, pdb_asym_id_2 = confidence_tools.compute_distances_from_pdb(
#             pdb_structure_2
#         )

#         interface_contact_mask_1 = custom_confidence.compute_interface_contact_mask(
#             distances=pdb_distances_1, asym_id=pdb_asym_id_1, dist=8
#         )
#         interface_contact_mask_2 = custom_confidence.compute_interface_contact_mask(
#             distances=pdb_distances_2, asym_id=pdb_asym_id_2, dist=8
#         )

#         plot_matrix(
#             os.path.join(
#                 args.outfolder, "predicted_tm_term_interface_contact_mask.png"
#             ),
#             predicted_tm_term_1 * interface_contact_mask_1,
#             predicted_tm_term_2 * interface_contact_mask_2,
#             limitA,
#             limitB,
#             name1=f"predicted_tm_term_{args.name1}",
#             name2=f"predicted_tm_term_{args.name2}",
#         )

#         interface_residues_mask_1 = custom_confidence.compute_interface_residues_mask(
#             distances=pdb_distances_1, asym_id=pdb_asym_id_1, dist=8
#         )
#         interface_residues_mask_2 = custom_confidence.compute_interface_residues_mask(
#             distances=pdb_distances_2, asym_id=pdb_asym_id_2, dist=8
#         )

#         plot_matrix(
#             os.path.join(
#                 args.outfolder, "predicted_tm_term_interface_residues_mask.png"
#             ),
#             predicted_tm_term_1 * interface_residues_mask_1,
#             predicted_tm_term_2 * interface_residues_mask_2,
#             limitA,
#             limitB,
#             name1=f"predicted_tm_term_{args.name1}",
#             name2=f"predicted_tm_term_{args.name2}",
#         )

#     logging.info("Plotting pae.")
#     pae1 = prediction_result_1["predicted_aligned_error"]
#     pae2 = prediction_result_2["predicted_aligned_error"]
#     plotfilepath = os.path.join(args.outfolder, "pae.png")
#     plot_matrix(
#         plotfilepath,
#         pae1,
#         pae2,
#         limitA,
#         limitB,
#         name1=f"pae_{args.name1}",
#         name2=f"pae_{args.name2}",
#     )

#     logging.info("Plotting plddt.")
#     plot_vec(
#         f"{args.outfolder}/plddt.png",
#         prediction_result_1["plddt"],
#         prediction_result_2["plddt"],
#         limitA=limitA,
#         name1=args.name1,
#         name2=args.name2,
#     )

#     logging.info("Plotting ptm alignment.")
#     tm_alignment_1 = custom_confidence.predicted_tm_alignment(
#         logits=prediction_result_1["predicted_aligned_error_raw"]["logits"],
#         breaks=prediction_result_1["predicted_aligned_error_raw"]["breaks"],
#     )
#     tm_alignment_2 = custom_confidence.predicted_tm_alignment(
#         logits=prediction_result_2["predicted_aligned_error_raw"]["logits"],
#         breaks=prediction_result_2["predicted_aligned_error_raw"]["breaks"],
#     )
#     plot_vec(
#         f"{args.outfolder}/tm_alignment.png",
#         tm_alignment_1,
#         tm_alignment_2,
#         limitA=limitA,
#         name1=args.name1,
#         name2=args.name2,
#     )

#     logging.info("Plotting iptm alignment.")
#     itm_alignment_1 = custom_confidence.predicted_tm_alignment(
#         logits=prediction_result_1["predicted_aligned_error_raw"]["logits"],
#         breaks=prediction_result_1["predicted_aligned_error_raw"]["breaks"],
#         asym_id=prediction_result_1["predicted_aligned_error_raw"]["asym_id"],
#         interface=True,
#     )
#     itm_alignment_2 = custom_confidence.predicted_tm_alignment(
#         logits=prediction_result_2["predicted_aligned_error_raw"]["logits"],
#         breaks=prediction_result_2["predicted_aligned_error_raw"]["breaks"],
#         asym_id=prediction_result_2["predicted_aligned_error_raw"]["asym_id"],
#         interface=True,
#     )
#     plot_vec(
#         f"{args.outfolder}/interface_tm_alignment.png",
#         itm_alignment_1,
#         itm_alignment_2,
#         limitA=limitA,
#         name1=args.name1,
#         name2=args.name2,
#     )

#     if os.path.exists(pdb_file_path_1) and os.path.exists(pdb_file_path_2):
#         logging.info("Calculating custom iptm scores.")
#         interface_contacts_1 = custom_confidence.compute_interface_contacts(
#             interface_contact_mask_1
#         )
#         interface_contacts_2 = custom_confidence.compute_interface_contacts(
#             interface_contact_mask_2
#         )

#         plot_vec(
#             f"{args.outfolder}/interface_contacts.png",
#             interface_contacts_1,
#             interface_contacts_2,
#             limitA=limitA,
#             name1=args.name1,
#             name2=args.name2,
#         )

#         iptm_interface_contacts_only_1 = np.asarray(
#             itm_alignment_1[(itm_alignment_1 * interface_contacts_1).argmax()]
#         )
#         iptm_interface_contacts_only_2 = np.asarray(
#             itm_alignment_2[(itm_alignment_2 * interface_contacts_2).argmax()]
#         )
#         metrics["iptm_interface_contacts_only_1"] = iptm_interface_contacts_only_1
#         metrics["iptm_interface_contacts_only_2"] = iptm_interface_contacts_only_2

#         plot_vec(
#             f"{args.outfolder}/interface_tm_alignment_restricted.png",
#             itm_alignment_1 * interface_contacts_1,
#             itm_alignment_2 * interface_contacts_2,
#             limitA=limitA,
#             name1=args.name1,
#             name2=args.name2,
#         )

#         mean_iptm_interface_contacts_only_1 = np.asarray(
#             np.mean(itm_alignment_1 * interface_contacts_1)
#         )
#         mean_iptm_interface_contacts_only_2 = np.asarray(
#             np.mean(itm_alignment_2 * interface_contacts_2)
#         )
#         metrics["mean_iptm_interface_contacts_only_1"] = (
#             mean_iptm_interface_contacts_only_1
#         )
#         metrics["mean_iptm_interface_contacts_only_2"] = (
#             mean_iptm_interface_contacts_only_2
#         )

#     # Compute MMalign

#     pdb_file_path_1 = os.path.join(
#         os.path.split(pickle_file_path_1)[0],
#         "_".join(
#             ["unrelaxed"] + os.path.split(pickle_file_path_1)[1].split("_")[1:]
#         ).split(".")[0]
#         + ".pdb",
#     )
#     pdb_file_path_2 = os.path.join(
#         os.path.split(pickle_file_path_2)[0],
#         "_".join(
#             ["unrelaxed"] + os.path.split(pickle_file_path_2)[1].split("_")[1:]
#         ).split(".")[0]
#         + ".pdb",
#     )

#     if (
#         os.path.exists(pdb_file_path_1)
#         and os.path.exists(pdb_file_path_2)
#         and args.MMalign_exe is not None
#     ):
#         logging.info("Calculation MMalign score.")
#         logging.info(f"Model 1: {os.path.split(pdb_file_path_1)[1]}")
#         logging.info(f"Model 2: {os.path.split(pdb_file_path_2)[1]}")

#         MMalign_result = MMalign_wrapper.run_MMalign(
#             args.MMalign_exe, pdb_file_path_1, pdb_file_path_2
#         )
#         #'PDBchain1', 'PDBchain2', 'TM1', 'TM2', 'RMSD', 'ID1', 'ID2', 'IDali', 'L1', 'L2', 'Lali'
#         metrics["TM1"] = MMalign_result["TM1"]
#         metrics["TM2"] = MMalign_result["TM2"]

#     metrics_filename = os.path.join(args.outfolder, "metrics.csv")
#     write_results(metrics_filename, metrics)

#     # Plot dgram and pae
#     logging.info("Plotting dgram and pae.")
#     dist_1, pae_1 = dgram2dmap_mod.load_results(pickle_file_path_1, interpolate=True)
#     dgram2dmap_mod.compare_to_native(
#         f"{args.outfolder}/comp_to_native_1.png",
#         pdb_file_path_1,
#         dist_1,
#         limitA,
#         limitB,
#     )
#     dgram2dmap_mod.plot_distances(
#         f"{args.outfolder}/dmap_1.png", dist_1, pae_1, limitA, limitB
#     )

#     dist_2, pae_2 = dgram2dmap_mod.load_results(pickle_file_path_2, interpolate=True)
#     dgram2dmap_mod.compare_to_native(
#         f"{args.outfolder}/comp_to_native_2.png",
#         pdb_file_path_2,
#         dist_2,
#         limitA,
#         limitB,
#     )
#     dgram2dmap_mod.plot_distances(
#         f"{args.outfolder}/dmap_2.png", dist_2, pae_2, limitA, limitB
#     )


# if __name__ == "__main__":
#     main()

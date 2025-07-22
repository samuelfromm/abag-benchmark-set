import confidence_tools
from typing import Optional
import pickle
import os
import sys
import argparse




def add_arguments(parser):
    parser.add_argument(
        "--results_pkl",
        help="path to results pkl file",
        type=str,
    )
    parser.add_argument(
        "--asym_id_pkl",
        help="path to pkl file containing the asymmetric ID. If not provided, it is assumed that the results pkl file containts the asymmetric unit ID.",
        type=str,
        default=None,
    )

    


def main():

    parser = argparse.ArgumentParser(
    description="Assemble new features for ligand receptor pair."
    )
    add_arguments(parser)
    args = parser.parse_args()

    results = confidence_tools.load_data_from_pkl(args.results_pkl)

    try:
        aligned_confidence_probs = results['aligned_confidence_probs']
    except:
        raise Exception("ERROR: 'aligned_confidence_probs' not present in results.")

    try:
        asym_id = results['asym_id']
    except:
        try:
            additional_results = confidence_tools.load_data_from_pkl(args.asym_id_pkl)
            asym_id = additional_results['asym_id']
        except:
            raise Exception("ERROR: 'asym_id' not present in results.")


    confidence_metrics = {}
    confidence_metrics['ptm'] = confidence_tools.calculate_predicted_tm_score(
        aligned_confidence_probs=aligned_confidence_probs,
        interface=False,
        residue_weights=None)
    confidence_metrics['iptm'] = confidence_tools.calculate_predicted_tm_score(
        aligned_confidence_probs=aligned_confidence_probs,
        asym_id=asym_id,
        interface=True,
        residue_weights=None)
    confidence_metrics['ranking_confidence'] = 0.8 * confidence_metrics['iptm'] + 0.2 * confidence_metrics['ptm']

    outstring = [f"{key}:{value}" for key,value in confidence_metrics.items()]
    outstring=','.join(outstring)
    print(outstring)



if __name__ == '__main__':
    main()


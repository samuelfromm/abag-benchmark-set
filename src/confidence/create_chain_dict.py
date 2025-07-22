import confidence_tools
from typing import Optional
import pickle
import os
import sys
import argparse




def add_arguments(parser):
    parser.add_argument(
        "-A",
        "--alignment",
        nargs=2,
        type=str,
        default=None,
        help='native chains and model chains superimposed onto each other. Chains are seperated by ":". Native chains come first. Example: "A:B:C" "L:H:A"',
    )

    


def main():

    parser = argparse.ArgumentParser(
    description="Create chain ID dictionary."
    )
    add_arguments(parser)
    args = parser.parse_args()

    native_chain_ids = args.alignment[0].split(":")
    model_chain_ids = args.alignment[1].split(":")
    if "" in native_chain_ids or "" in model_chain_ids:
        # example: A:B::C E:F:G:: -> use alignment: ABC EFG
        #print(
        #    "WARNING: Not all chains are aligned. Replacing alignment with basic alignment order."
        #)
        native_chain_ids = [chain_id for chain_id in native_chain_ids if chain_id != ""]
        model_chain_ids = [chain_id for chain_id in model_chain_ids if chain_id != ""]
    assert len(native_chain_ids) == len(model_chain_ids)

    model_chain_id_dict = {model_chain:idx for idx,model_chain in enumerate(sorted(model_chain_ids),start=1)}

    conversion_dict = {native_chain_ids[i]: model_chain_ids[i] for i in range(len(native_chain_ids))}
    conversion_dict = dict(sorted(conversion_dict.items()))

    outstring=''
    outstring+='model_chain_id_dict:'+' '.join([f"{key}:{value}" for key,value in model_chain_id_dict.items])+';'
    outstring+='conversion_dict:'+' '.join([f"{key}:{value}" key,value in model_chain_id_dict.items])+';'

    print(outstring)



if __name__ == '__main__':
    main()


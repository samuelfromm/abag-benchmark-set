import pickle
import json
import sys
import numpy as np


def read_pkl(pkl_file):
    # Load data from the pickle file
    with open(pkl_file, "rb") as f:
        data = pickle.load(f)

    print(data.keys())


# Check for input arguments
if len(sys.argv) != 2:
    print("Usage: python read_pkl.py input.pkl")
else:
    pkl_file = sys.argv[1]
    read_pkl(pkl_file)

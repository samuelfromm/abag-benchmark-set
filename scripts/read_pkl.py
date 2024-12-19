import pickle
import sys

def read_pkl(pkl_file, keys=None):
    # Load data from the pickle file
    with open(pkl_file, "rb") as f:
        data = pickle.load(f)

    print("Available keys:", data.keys())

    if keys:
        for key in keys:
            if key in data:
                print(f"Data for key '{key}': {data[key]}")
            else:
                print(f"Key '{key}' not found in the data.")


# Check for input arguments
if len(sys.argv) < 2:
    print("Usage: python read_pkl.py input.pkl [key1 key2 ...]")
else:
    pkl_file = sys.argv[1]
    keys = sys.argv[2:] if len(sys.argv) > 2 else None
    read_pkl(pkl_file, keys)
import pickle
import sys
import json
import os
import numpy as np

def read_file(file_path, keys=None):
    # Check the file extension to determine the format
    _, file_extension = os.path.splitext(file_path)
    
    if file_extension == ".pkl":
        # Load data from the pickle file
        with open(file_path, "rb") as f:
            data = pickle.load(f)
    elif file_extension == ".json":
        # Load data from the JSON file
        with open(file_path, "r") as f:
            data = json.load(f)
    else:
        print(f"Unsupported file format: {file_extension}. Please use a .pkl or .json file.")
        return
    
    # Print available keys in the data
    if isinstance(data, dict):
        print("Available keys:", data.keys())
    else:
        print("The data is not a dictionary. It might not have keys to inspect.")
        return

    # Print data for specified keys
    if keys:
        for key in keys:
            if key in data:
                print(f"Data for key '{key}': {data[key]} (type: {type(data[key])})")
                try: 
                    print(f"Shape: {data[key].shape}")
                except:
                    pass
            else:
                print(f"Key '{key}' not found in the data.")

# Check for input arguments
if len(sys.argv) < 2:
    print("Usage: python read_file.py input_file.[pkl|json] [key1 key2 ...]")
else:
    file_path = sys.argv[1]
    keys = sys.argv[2:] if len(sys.argv) > 2 else None
    read_file(file_path, keys)
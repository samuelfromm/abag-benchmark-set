import pickle
import json
import sys

def pkl_to_json(pkl_file, json_file):
    # Load data from the pickle file
    with open(pkl_file, 'rb') as f:
        data = pickle.load(f)
    
    # Save data to a JSON file
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=4)
    
    print(f"Data successfully converted from {pkl_file} to {json_file}")

# Check for input arguments
if len(sys.argv) != 3:
    print("Usage: python pkl_to_json.py input.pkl output.json")
else:
    pkl_file = sys.argv[1]
    json_file = sys.argv[2]
    pkl_to_json(pkl_file, json_file)
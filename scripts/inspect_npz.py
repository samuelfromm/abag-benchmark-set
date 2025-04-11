import sys
import numpy as np

def inspect_npz(file_path):
    try:
        # Load the .npz file
        data = np.load(file_path)
        
        print(f"Inspecting contents of: {file_path}\n")
        print(f"File contains {len(data.files)} arrays:\n")

        # Iterate over each item in the .npz file
        for name in data.files:
            array = data[name]
            print(f"Array name: {name}")
            print(f"  Shape: {array.shape}")
            print(f"  Dtype: {array.dtype}")
            print()

    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred while inspecting the file: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python inspect_npz.py <path_to_npz_file>")
    else:
        file_path = sys.argv[1]
        inspect_npz(file_path)
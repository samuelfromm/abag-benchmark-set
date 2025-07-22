import torch
import openfold.data.data_transforms
import openfold.np.protein
from openfold.utils.rigid_utils import Rigid
from typing import Optional


def read_pdb(pdb_path: str):
    """Read a pdb file and return the content as a string."""
    with open(pdb_path, "r") as pdb_file:
        pdb_str = pdb_file.read()
    return pdb_str


def protein_from_pdb(pdb_path: str):

    protein = openfold.np.protein.from_pdb_string(read_pdb(pdb_path))

    return protein

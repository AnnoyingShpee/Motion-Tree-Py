import gemmi
import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import cdist
"""
Gemmi follows a hierarchy:
Structure -> Model -> Chain -> [ResidueSpan] -> Residue -> Atom
"""


class Protein:
    def __init__(self, name, file_path: str, chain: str = "A"):
        self.name = name
        self.file_path = file_path
        self.chain_param: str = chain  # The chain specified in parameter input for the program
        self.utilised_res_nums = None
        # Dataframe that stores the coordinates of the utilised atoms of the residues. Only Ca atoms or backbone atoms.
        self.utilised_atoms_coords = None
        self.distance_matrix = None

    def get_distance_matrix(self):
        self.distance_matrix = cdist(self.utilised_atoms_coords, self.utilised_atoms_coords, metric="euclidean")

    def get_structure(self):
        return gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)

    def get_model(self):
        structure = gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)
        return structure[0]

    def get_chain(self):
        # There is usually only one model in the structure
        structure = gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)
        return structure[0][self.chain_param]

    def get_polymer(self):
        structure = gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)
        return structure[0][self.chain_param].get_polymer()

    def print_chain(self):
        print(f"{self.get_structure().name}({self.chain_param}) - {self.utilised_atoms_coords.shape}")
        print(f"{self.get_structure().name}({self.chain_param}) - \n{self.utilised_atoms_coords}")


import gemmi
import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances
"""
Gemmi follows a hierarchy:
Structure -> Model -> Chain -> [ResidueSpan] -> Residue -> Atom
"""


class Protein:
    def __init__(self, name, file_path: str, chain: str = "A", atom_type: str = "ca"):
        self.name = name
        self.file_path = file_path
        self.chain_param: str = chain  # The chain specified in parameter input for the program
        self.atom_type: str = atom_type
        self.utilised_res_nums = []
        # Dataframe that stores the coordinates of the utilised atoms of the residues. Only Ca atoms or backbone atoms.
        self.main_atoms_coords = None
        if atom_type == "ca":
            atoms_coords, self.utilised_res_nums = self.remove_no_ca_residues()
            self.main_atoms_coords = np.asarray(atoms_coords)
        elif atom_type == "backbone":
            atoms_coords, self.utilised_res_nums = self.remove_no_backbone_residues()
            self.main_atoms_coords = np.asarray(atoms_coords)
        else:
            print("No atom type given. Default to Ca atoms.")
            atoms_coords, self.utilised_res_nums = self.remove_no_ca_residues()
            self.main_atoms_coords = np.asarray(atoms_coords)
        self.distance_matrix = None
        self.get_distance_matrix()

    def remove_no_ca_residues(self):
        """
        Checks if each residue has a Ca atom. If no Ca atom present in residue, remove residue.
        :return:
        """
        utilised_residue_nums = []
        atom_coords = []
        chain = self.get_chain()
        for r in range(len(chain)):
            for a in chain[r]:
                if a.name == "CA":
                    atom_coord: gemmi.Atom = a
                    utilised_residue_nums.append(r)
                    atom_coords.append(atom_coord.pos.tolist())
                    break

        return atom_coords, utilised_residue_nums

    def remove_no_backbone_residues(self):
        utilised_residue_nums = []
        atom_coords = []
        chain = self.get_chain()
        for r in range(len(chain)):
            for a in chain[r]:
                if a.name == "N" or a.name == "CA" or a.name == "C":
                    atom_coord: gemmi.Atom = a
                    if r not in utilised_residue_nums:
                        utilised_residue_nums.append(r)
                    atom_coords.append(atom_coord.pos.tolist())

        return atom_coords, utilised_residue_nums

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

    def get_distance_matrix(self):
        self.distance_matrix = pairwise_distances(self.main_atoms_coords, self.main_atoms_coords)

    def print_chain(self):
        print(f"{self.get_structure().name}({self.chain_param}) - {self.main_atoms_coords.shape}")
        print(f"{self.get_structure().name}({self.chain_param}) - \n{self.main_atoms_coords}")


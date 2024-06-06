import gemmi
import pandas as pd
from scipy.spatial import distance_matrix
"""
Gemmi follows a hierarchy:
Structure -> Model -> Chain -> Residue -> Atom
"""

backbone_atoms = ["N", "CA", "C"]


class Protein:
    def __init__(self, file_path: str, chain: str = "A", atom_type: str = "ca"):
        self.id: str = self.structure.name  # ID of the protein structure
        self.file_path = file_path
        self.chain_param: str = chain  # The chain specified in parameter input for the program
        self.atom_type: str = atom_type
        self.structure: gemmi.Structure = gemmi.read_structure(file_path, format=gemmi.CoorFormat.Pdb)
        # There is usually only one model in the structure
        self.chain: gemmi.Chain = self.structure[0][self.chain_param]
        self.residue_span: gemmi.ResidueSpan = self.chain.get_polymer()
        self.chain_atoms = None
        self.slide_window_residues_indices = None
        self.slide_window_residues = []
        self.distance_matrix = None
        self.get_backbone_atoms()

    def get_backbone_atoms(self):
        """
        Gets the backbone atoms of each residue [(N, CA, C) or (CA only)]
        :return: A 2D array of residue atoms
        """
        atoms = []
        if self.atom_type == "backbone":
            for res in self.residue_span:
                atoms.append([res.sole_atom(a) for a in backbone_atoms])
        elif self.atom_type == "ca":
            for res in self.residue_span:
                atoms.append(res.sole_atom("CA"))
        self.chain_atoms = np.array(atoms)

    def recreate_structure(self):
        """
        To reconstruct the protein to its original structure after any changes to it.
        :return:
        """
        self.slide_window_residues = []
        self.structure: gemmi.Structure = gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)
        self.chain: gemmi.Chain = self.structure[0][self.chain_param]
        self.residue_span: gemmi.ResidueSpan = self.chain.get_polymer()
        for i in range(self.slide_window_residues_indices[0], self.slide_window_residues_indices[1]):
            self.slide_window_residues.append(self.chain[i])
        self.id = self.structure.name

    def get_structure(self):
        return gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)

    def get_model(self):
        structure = gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)
        return structure[0]

    def get_chain(self):
        structure = gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)
        return structure[0][self.chain_param]

    def get_polymer(self):
        structure = gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)
        return structure[0][self.chain_param].get_polymer()

    def get_slide_window_result(self):
        structure = gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)
        slide_window = structure[0][self.chain_param].get_polymer()
        for i in range(len(slide_window) - self.slide_window_residues_indices[1]):
            slide_window.__delitem__(-1)
        for i in range(self.slide_window_residues_indices[0]):
            slide_window.__delitem__(0)
        return slide_window

    def get_distance_matrix(self):

        pass

    def print_chain(self):
        print(f"{self.id}({self.chain_param}) - {self.chain_atoms.shape}")
        print(f"{self.id}({self.chain_param}) - {self.chain_atoms}")

    def print_slide_window_residues(self, n=None):
        if n is None or n > len(self.slide_window_residues):
            n = len(self.slide_window_residues)
        print(f"slide_window_residue_1 shape = {len(self.slide_window_residues)}")
        print(f"slide_window_residues_1[0:{n}] = ")
        print(self.slide_window_residues)
        # [print(f"{self.slide_window_residues[r]}") for r in range(0, n)]



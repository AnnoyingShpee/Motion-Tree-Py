import gemmi
import numpy as np
from scipy.spatial.distance import cdist
"""
Gemmi follows a hierarchy:
Structure -> Model -> Chain -> [ResidueSpan] -> Residue -> Atom
"""


class Protein:
    def __init__(self, input_path, code, chain: str = "A", is_dyndom=False):
        self.input_path = input_path
        self.code = code
        self.file_path = f"{input_path}/{code}.pdb"
        self.chain_param: str = chain  # The chain specified in parameter input for the program
        self.is_dyndom = is_dyndom
        # The indices of the residues that are utilised in the protein chain
        self.utilised_res_indices = None
        # Dataframe that stores the coordinates of the utilised atoms of the residues. Only Ca atoms or backbone atoms.
        self.utilised_atoms_coords = None
        self.distance_matrix = None
        structure = gemmi.read_pdb(self.file_path)

    def get_distance_matrix(self):
        # print(self.utilised_atoms_coords.shape)
        self.distance_matrix = cdist(self.utilised_atoms_coords, self.utilised_atoms_coords, metric="euclidean")
        # print(self.distance_matrix.shape)

    def get_structure(self):
        return gemmi.read_pdb(self.file_path)

    def get_model(self):
        structure = gemmi.read_pdb(self.file_path)
        if self.is_dyndom:
            return structure[1]
        return structure[0]

    def get_chain(self):
        # There is usually only one model in the structure
        structure = gemmi.read_pdb(self.file_path)
        if self.is_dyndom and self.chain_param == "B":
            return structure[1][self.chain_param]
        else:
            return structure[0][self.chain_param]

    def get_polymer(self):
        structure = gemmi.read_pdb(self.file_path)
        return structure[0][self.chain_param].get_polymer()

    def get_polymer_entity(self):
        structure: gemmi.Structure = gemmi.read_pdb(self.file_path)
        polymer = self.get_polymer()
        return structure.get_entity_of(polymer)

    def get_residue_nums(self, indices, utilised=True):
        if self.is_dyndom:
            polymer = self.get_chain()
        else:
            polymer = self.get_polymer()
        if utilised:
            return [polymer[i].seqid.num for i in self.utilised_res_indices[indices]]
        else:
            temp = np.delete(self.utilised_res_indices, indices)
            return [polymer[i].seqid.num for i in temp]

    def print_chain(self):
        print(f"{self.get_structure().name}({self.chain_param}) - {self.utilised_atoms_coords.shape}")
        print(f"{self.get_structure().name}({self.chain_param}) - \n{self.utilised_atoms_coords}")

    def print_dist_mat(self):
        row = [str(self.utilised_res_indices[i]).ljust(4) for i in range(self.distance_matrix.shape[1])]
        print("[".ljust(3), " ", " ".join(row))
        for i in range(self.distance_matrix.shape[0]):
            row = [str(round(j, 1)).ljust(4) for j in self.distance_matrix[i]]
            print(str(i).ljust(3), " ", " ".join(row))
        print("]")


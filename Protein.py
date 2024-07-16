import gemmi
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
        # The indices of the residues that are utilised in the protein chain
        self.utilised_res_indices = None
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

    def get_residue_nums(self, indices):
        polymer = self.get_polymer()
        return [polymer[i].seqid.num for i in self.utilised_res_indices[indices]]

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


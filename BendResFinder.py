import gemmi
import Protein
import numpy as np


class BendResFinder:
    def __init__(self, protein_1, protein_2, nodes):
        self.protein_1: Protein = protein_1
        self.protein_2: Protein = protein_2
        self.protein_1_slide_indices = None
        self.protein_2_slide_indices = None
        self.nodes: dict = nodes
        self.window_size = 5
        self.fitted_protein_2 = None
        self.chain_superimpose_result = None
        self.slide_window_superimpose_results = None
        # Bending residues indices
        self.bending_residues = {}

    def superimpose_chains(self):
        """
        Superimposes the entire chain of Protein 2 onto Protein 1 using the backbone atoms to get Protein 2 in the same
        coordinate space of Protein 1. Saves the Protein 2 chain after transformation into fitting_protein_2.
        """
        ptype = self.protein_1.get_polymer().check_polymer_type()

        self.chain_superimpose_result: gemmi.SupResult = gemmi.calculate_superposition(self.protein_1.get_polymer(),
                                                                                       self.protein_2.get_polymer(),
                                                                                       ptype, gemmi.SupSelect.MainChain)
        self.fitted_protein_2 = self.protein_2.get_chain()
        polymer = self.fitted_protein_2.get_polymer()
        polymer.transform_pos_and_adp(self.chain_superimpose_result.transform)

    def get_slide_window_start_end_indices(self):
        """
        Gets the start and end indices of the sliding windows using the index of the middle residue in each window.
        If window size is 5 and the length of the utilised residues is 431 (Index 430), the start and end indices are
        2 and 428 respectively.
        :return:
        """
        window_size = self.window_size
        # The index of the middle residue in the sliding window
        window_mid_index = start_index = (window_size - 1) // 2
        final_index = len(self.protein_1.utilised_res_indices) + 1 - (window_size - window_mid_index)
        self.protein_1_slide_indices = (start_index, final_index)
        self.protein_2_slide_indices = (start_index, final_index)

    def sliding_window_superimpose_residues(self):
        """
        Slides a window of a specified size over both protein chains. The backbone atoms of the residues in the
        sliding windows are superimposed to get the superposition results of the middle residue of each sliding window.
        :return:
        """
        self.slide_window_superimpose_results = []
        fitting_protein_polymer = self.protein_1.get_polymer()  # The protein that will superimpose onto Protein S.
        target_protein_polymer = self.fitted_protein_2.get_polymer()
        # For each window, get the residues' backbone atoms and get the superposition results.
        for r in range(len(self.protein_1.utilised_residues_indices) - self.window_size + 1):  # Number of iterations of the window
            # Initialise an empty list that stores the backbone atoms' coordinates of the target protein
            # (The fitted protein 2)
            target_protein_polymer_atoms_pos = []
            # Initialise an empty list that stores the backbone atoms' coordinates of the fitting protein
            # (The protein 1 that will fit onto the fitted protein 2)
            fitting_protein_polymer_atoms_pos = []
            # For each residue in the window,
            for i in range(r, r+self.window_size):
                fitting_index = self.protein_1.utilised_res_indices[i]
                target_index = self.protein_2.utilised_res_indices[i]
                for a in fitting_protein_polymer[fitting_index]:
                    if a.name == "CA":
                        fitting_protein_polymer_atoms_pos.append(a.pos)
                        break
                for a in target_protein_polymer[target_index]:
                    if a.name == "CA":
                        target_protein_polymer_atoms_pos.append(a.pos)
                        break
            # Superimpose the atoms and append the result to list
            self.slide_window_superimpose_results.append(gemmi.superpose_positions(target_protein_polymer_atoms_pos,
                                                                                   fitting_protein_polymer_atoms_pos))

    def superimpose_chains(self):
        self.fitted_protein_2 = self.protein_2.get_chain()
        fitted_protein_2_polymer = self.fitted_protein_2.get_polymer()
        print(self.fitted_protein_2[50]['CA'].pos)
        ptype = fitted_protein_2_polymer.check_polymer_type()

        self.chain_superimpose_result: gemmi.SupResult = gemmi.calculate_superposition(self.protein_1.get_polymer(),
                                                                                       fitted_protein_2_polymer,
                                                                                       ptype, gemmi.SupSelect.CaP)
        fitted_protein_2_polymer.transform_pos_and_adp(self.chain_superimpose_result.transform)
        print(self.fitted_protein_2[50]['CA'].pos)

    def slide_window_superimpose_residues(self):
        self.slide_window_superimpose_results = []


    def get_rotation_vectors(self):

        return

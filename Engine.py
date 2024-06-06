import numpy as np
import pandas as pd
import gemmi
from Protein import Protein
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.transform import Rotation


class Engine:
    def __init__(self, first_pdb: str, second_pdb: str, commands=None, params=None):
        # Default parameters to be used if not given input
        if commands is None:
            self.commands = {
                "chain1id": "A",
                "chain2id": "A",
                "atoms": "ca"
            }
        else:
            self.commands = commands
        self.parameters = params
        # Initialise proteins
        self.protein_1: Protein = Protein(first_pdb, self.commands["chain1id"], self.commands["atoms"])
        self.protein_2: Protein = Protein(second_pdb, self.commands["chain2id"], self.commands["atoms"])
        self.protein_1_atoms = None
        self.protein_2_atoms = None
        # print(len(self.protein_1.get_polymer()))
        # print(len(self.protein_2.get_polymer()))
        self.aligned_protein_1 = None
        self.main_atoms = []
        if self.commands["atoms"] == "backbone":
            self.main_atoms = ["N", "CA", "C"]
        elif self.commands["atoms"] == "ca":
            self.main_atoms = ["CA"]
        # A object containing superimposition information between the proteins
        self.chain_superimpose_result = None
        # Bending residues indices
        self.bending_residues = set()

    def run(self):
        """
        Runs the entire program. To make it simple, all operations happen based on Protein 1 conforming to Protein 2.
        :return:
        """
        running = True
        while running:
            # Superimposes Protein 2 onto Protein 1, so they are in the same "space". The superimposed protein will be
            # called Protein S.
            self.superimpose_chains()

            agglo_cluster = AgglomerativeClustering(linkage=self.parameters["linkage"])

            running = False
        return True

    def superimpose_chains(self):
        """
        Superimposes the entire chain of Protein 1 unto Protein 2 using the backbone atoms.
        :return: Protein 2 chain after transformation
        """
        ptype = self.protein_1.residue_span.check_polymer_type()
        atom_selection = None
        if self.commands["atoms"] == "ca":
            atom_selection = gemmi.SupSelect.CaP
        elif self.commands["atoms"] == "backbone":
            atom_selection = gemmi.SupSelect.MainChain
        else:
            atom_selection = gemmi.SupSelect.All
        fit_1_to_2_result: gemmi.SupResult = gemmi.calculate_superposition(self.protein_2.get_polymer(),
                                                                           self.protein_1.get_polymer(),
                                                                           ptype, atom_selection)

        # print(f"Count = {fit_1_to_2_result.count}")
        # print(f"Center 1 = {fit_1_to_2_result.center1}")
        # print(f"Center 2 = {fit_1_to_2_result.center2}")
        # print(f"RMSD = {fit_1_to_2_result.rmsd}")
        # print(f"Rotation = {fit_1_to_2_result.transform.mat}")
        # print(f"Translation = {fit_1_to_2_result.transform.vec}")
        self.chain_superimpose_result = fit_1_to_2_result
        self.aligned_protein_1: gemmi.ResidueSpan = self.protein_1.get_polymer()
        self.aligned_protein_1.transform_pos_and_adp(fit_1_to_2_result.transform)

    def calculate_distance_between_points(self):
        """
        Calculates the Euclidean distance between 2 atoms.
        :return:
        """
        return

    def calculate_distance_difference_matrix(self):
        """
        Uses the Euclidean distances from the 2 chains to calculate the distance difference matrix.
        :return:
        """
        return

    def get_utilised_residues(self, indices):
        """
        Gets all residues that were used in the middle of each sliding window
        :param indices:
        :return:
        """
        chain_1 = self.protein_1.get_polymer()
        chain_2 = self.protein_2.get_polymer()
        for i in range(indices[0], indices[1]):
            self.protein_1.slide_window_residues.append(chain_1[i])
            self.protein_2.slide_window_residues.append(chain_2[i])

    def determine_bending_residues(self):
        fixed_domain = self.clusterer.domains[self.clusterer.fixed_domain]
        fixed_domain_segments = fixed_domain.segments
        fixed_domain_rot_vecs = self.rotation_vecs[fixed_domain_segments[0][0]:fixed_domain_segments[0][1]]

        # Get the rotation vectors of the fixed domain
        for i in range(1, fixed_domain_segments.shape[0]):
            rot_vecs = self.rotation_vecs[fixed_domain_segments[i][0]:fixed_domain_segments[i][1]]
            fixed_domain_rot_vecs = np.append(fixed_domain_rot_vecs, rot_vecs, axis=0)
        # Calculate mean and standard deviation of the fixed domain rotation vectors
        fixed_domain_mean = np.mean(fixed_domain_rot_vecs, axis=0)
        fixed_domain_std = np.std(fixed_domain_rot_vecs)
        fixed_domain_centered_vecs = fixed_domain_rot_vecs - fixed_domain_mean
        fixed_domain_covar = np.cov(fixed_domain_centered_vecs.T)
        fixed_domain_inv_covar = np.linalg.inv(fixed_domain_covar)
        fixed_domain_var = np.diag(fixed_domain_covar)
        # print(f"Fixed Domain Mean = {fixed_domain_mean}")
        # print(f"Fixed Domain STD = {fixed_domain_std}")
        # print(f"Fixed Domain Var = {fixed_domain_var}")
        # print(f"Fixed Domain Covariance = \n{fixed_domain_covar}")
        # print(f"Fixed Domain Inverse Covariance = \n{fixed_domain_inv_covar}")

        # For each dynamic domain,
        for domain in self.clusterer.domains:
            if domain.domain_id == self.clusterer.fixed_domain:
                continue
            print(f"Domain {domain.domain_id}")
            dyn_dom_segments = domain.segments
            dyn_dom_rot_vecs = self.rotation_vecs[dyn_dom_segments[0][0]:dyn_dom_segments[0][1]]

            for i in range(1, dyn_dom_segments.shape[0]):
                rot_vecs = self.rotation_vecs[dyn_dom_segments[i][0]:dyn_dom_segments[i][1]]
                dyn_dom_rot_vecs = np.append(dyn_dom_rot_vecs, rot_vecs, axis=0)

            dyn_dom_mean = np.mean(dyn_dom_rot_vecs, axis=0)
            dyn_dom_std = np.std(dyn_dom_rot_vecs)
            dyn_dom_centered_vecs = dyn_dom_rot_vecs - dyn_dom_mean
            dyn_dom_covar = np.cov(dyn_dom_centered_vecs.T)
            dyn_dom_inv_covar = np.linalg.inv(dyn_dom_covar)
            dyn_dom_var = np.diag(dyn_dom_covar)

            # print(f"Dyn Domain Mean = {dyn_dom_mean}")
            # print(f"Dyn Domain STD = {dyn_dom_std}")
            # print(f"Dyn Domain Var = {dyn_dom_var}")
            # print(f"Dyn Domain Covariance = \n{dyn_dom_covar}")
            # print(f"Dyn Domain Inverse Covariance = \n{dyn_dom_inv_covar}")

            # Calculate the indices of the previous and next residues for each segment of the dynamic domain.
            dyn_dom_prev_indices = dyn_dom_segments[:, 0] - 1
            dyn_dom_next_indices = dyn_dom_segments[:, 1] + 1

            # print("Fixed Domain segments: ")
            # print(fixed_domain.segments)
            # print("Dyn Dom segments:")
            # print(dyn_dom_segments)
            # print("Dyn Dom prev indices: ")
            # print(dyn_dom_prev_indices)
            # print("Dyn Dom next indices: ")
            # print(dyn_dom_next_indices)

            # Get the indices of the fixed domain segments that connects the fixed domain to the dynamic domain.
            # 1D Array of booleans where True means next index after fixed domain segment is dyn dom segment.
            fixed_next_is_dyn = np.in1d(fixed_domain.segments[:, 1], dyn_dom_prev_indices)
            fixed_next_is_dyn_ind = np.where(fixed_next_is_dyn)[0]
            # 1D Array of booleans where True means previous index before fixed domain segment is dyn dom segment.
            fixed_prev_is_dyn = np.in1d(fixed_domain.segments[:, 0], dyn_dom_next_indices)
            fixed_prev_is_dyn_ind = np.where(fixed_prev_is_dyn)[0]

            # Get the indices of the dynamic domain segments that connects the dynamic domain to the fixed domain.
            # 1D Array of booleans where True means next index after dyn dom segment is fixed domain segment.
            dyn_next_is_fixed = np.in1d(dyn_dom_next_indices, fixed_domain.segments[:, 0])
            dyn_next_is_fixed_ind = np.where(dyn_next_is_fixed)[0]
            # 1D Array of booleans where True means previous index before dyn dom segment is fixed domain segment.
            dyn_prev_is_fixed = np.in1d(dyn_dom_prev_indices, fixed_domain.segments[:, 1])
            dyn_prev_is_fixed_ind = np.where(dyn_prev_is_fixed)[0]

            print("Fixed domain segment next index is dyn dom segment: ")
            print(fixed_next_is_dyn)
            print(fixed_next_is_dyn_ind)
            print("Fixed domain segment prev index is dyn dom segment: ")
            print(fixed_prev_is_dyn)
            print(fixed_prev_is_dyn_ind)

            print("Dyn dom segment next index is fixed domain segment: ")
            print(dyn_next_is_fixed)
            print(dyn_next_is_fixed_ind)
            print("Dyn dom segment prev index is fixed domain segment: ")
            print(dyn_prev_is_fixed)
            print(dyn_prev_is_fixed_ind)

            p = 4.6

            # Go backwards through the fixed domain residues of the segments
            for segment_ind in fixed_next_is_dyn_ind:
                segment = fixed_domain.segments[segment_ind]
                for i in range(segment[1], segment[1] - 1, -1):
                    centered_vec = self.rotation_vecs[i] - fixed_domain_mean
                    q_value = centered_vec @ fixed_domain_inv_covar @ centered_vec
                    # print("Backward Fixed Q Value =", q_value)
                    if q_value > p:
                        self.bending_residues.add(i)
                    else:
                        break

            # Go forwards through the dyn dom residues of the segments
            for segment_ind in dyn_prev_is_fixed_ind:
                segment = domain.segments[segment_ind]
                for i in range(segment[1], segment[0] + 1):
                    centered_vec = self.rotation_vecs[i] - dyn_dom_mean
                    q_value = centered_vec @ dyn_dom_inv_covar @ centered_vec
                    # print("Forward Dyn Q Value =", q_value)
                    if q_value > p:
                        self.bending_residues.add(i)
                    else:
                        break

            # Go forwards through the fixed domain residues of the segments
            for segment_ind in fixed_prev_is_dyn_ind:
                segment = fixed_domain.segments[segment_ind]
                for i in range(segment[0], segment[1] + 1):
                    centered_vec = self.rotation_vecs[i] - fixed_domain_mean
                    q_value = centered_vec @ fixed_domain_inv_covar @ centered_vec
                    # print("Forward Fixed Q Value =", q_value)
                    if q_value > p:
                        self.bending_residues.add(i)
                    else:
                        break

            # Go backwards through the dyn dom residues of the segments
            for segment_ind in dyn_next_is_fixed_ind:
                segment = domain.segments[segment_ind]
                for i in range(segment[1], segment[0] - 1, -1):
                    centered_vec = self.rotation_vecs[i] - dyn_dom_mean
                    q_value = centered_vec @ dyn_dom_inv_covar @ centered_vec
                    # print("Backward Dyn Q Value =", q_value)
                    # if q_value < p or q_value > 1-p:
                    if q_value > p:
                        self.bending_residues.add(i)
                    else:
                        break

            print(self.bending_residues)
            print(len(self.bending_residues))

    def get_fixed_domain_transformations(self):
        slide_window_1 = self.protein_1.get_slide_window_result()
        slide_window_2 = self.protein_2.get_slide_window_result()
        coords_1 = []
        coords_2 = []
        fixed_domain_id = self.clusterer.fixed_domain
        for s in self.clusterer.domains[fixed_domain_id].segments:
            for i in range(s[0], s[1] + 1):
                for a in self.main_atoms:
                    coords_1.append(slide_window_1[i].sole_atom(a).pos)
                    coords_2.append(slide_window_2[i].sole_atom(a).pos)

        r: gemmi.SupResult = gemmi.superpose_positions(coords_1, coords_2)
        return r

    def determine_external_component_motion(self):

        return

    def print_chains_superimposed_result(self):
        # A gemmi.SupResult object containing superimposition information between 2 chains
        print(f"RMSD =                  {self.chain_superimpose_result.rmsd}")
        print(f"Count =                 {self.chain_superimpose_result.count}")
        print(f"Center 1 =              {self.chain_superimpose_result.center1}")
        print(f"Center 2 =              {self.chain_superimpose_result.center2}")
        print(f"Translation Vector =    {self.chain_superimpose_result.transform.vec}")
        print(f"Rotation Matrix =       {self.chain_superimpose_result.transform.mat}")

    def print_slide_window_superimpose_results(self, n=None):
        if n is None or n > len(self.slide_window_superimpose_results):
            n = len(self.slide_window_superimpose_results)
        print(f"slide_window_superimpose_result size = {len(self.slide_window_superimpose_results)}")
        for i in range(n):
            item: gemmi.SupResult = self.slide_window_superimpose_results[i]
            print(f"RMSD =                  {item.rmsd}")
            print(f"Count =                 {item.count}")
            print(f"Center 1 =              {item.center1}")
            print(f"Center 2 =              {item.center2}")
            print(f"Translation Vector =    {item.transform.vec}")
            print(f"Rotation Matrix =       {item.transform.mat}")

    def print_slide_window_residue_indices(self):
        print(f"slide_window_residue_indices = {self.protein_1.slide_window_residues_indices}")

    def print_unit_vectors(self, n=None):
        print(f"rotation_vecs shape = {self.unit_vectors.shape}")
        if n is None or n > self.unit_vectors.shape[0]:
            n = self.unit_vectors.shape[0]
        for i in range(n):
            print(
                f"{self.protein_1.slide_window_residues[i]} - [{self.unit_vectors[i][0]}, {self.unit_vectors[i][1]}, {self.unit_vectors[i][2]}]")


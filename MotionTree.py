import numpy as np
import gemmi
from copy import deepcopy
from timeit import default_timer
from statistics import mean
from Protein import Protein
from DataMngr import conn, ftp_files_to_disk, save_results_to_disk, write_info_file, write_to_pdb, \
    write_domains_to_pml  #, group_continuous_num, check_protein_exists, insert_protein_dist_mat, get_protein_dist_mat,
from difflib import SequenceMatcher


class MotionTree:
    def __init__(self, input_path, output_path, protein_1_name, chain_1, protein_2_name, chain_2,
                 spat_prox=7.0, small_node=5, clust_size=30, magnitude=5):
        self.input_path = input_path
        self.output_path = output_path
        self.protein_1_name = protein_1_name
        self.chain_1 = chain_1
        self.protein_2_name = protein_2_name
        self.chain_2 = chain_2
        self.spat_prox = spat_prox
        self.small_node = small_node
        self.clust_size = clust_size
        self.magnitude = magnitude
        # Initialise proteins
        self.protein_1 = None
        self.protein_2 = None
        self.res_nums = []
        self.similarity = 0
        # The starting distance difference matrix
        self.diff_dist_mat_init = None
        self.num_residues = None
        # Dictionary storing the clusters containing the index of the atoms
        # self.clusters = {i: [i] for i in range(self.diff_dist_mat_init.shape[0])}
        self.clusters = {}
        # The linkage matrix for building the dendrogram
        # self.link_mat = np.empty((self.diff_dist_mat_init.shape[0] - 1, 4))
        self.link_mat = None
        self.nodes = {}
        self.is_db_connected = True

        ftp_files_to_disk(self.input_path, self.protein_1_name, self.protein_2_name)

    def init_protein(self, protein_num):
        print("Init Protein", protein_num)
        if protein_num == 1:
            self.protein_1 = Protein(
                input_path=self.input_path,
                code=self.protein_1_name,
                chain=self.chain_1
            )
            return f"{self.protein_1_name} Initialised"
        else:
            self.protein_2 = Protein(
                input_path=self.input_path,
                code=self.protein_2_name,
                chain=self.chain_2
            )
            return f"{self.protein_2_name} Initialised"

    def check_sequence_identity(self):
        print("Checking Sequence")
        protein_1_polymer = self.protein_1.get_polymer()
        protein_2_polymer = self.protein_2.get_polymer()

        coords_1, coords_2, utilised_res_ind_1, utilised_res_ind_2, res_names_1, res_names_2 = self.get_ca_atoms_coords()
        # print(res_names_1)
        # print(res_names_2)
        correct_hit = 0
        # self.similarity = SequenceMatcher(None, res_names_1, res_names_2).ratio()
        for i in range(len(res_names_1)):
            if res_names_1[i] == res_names_2[i]:
                correct_hit += 1
        self.similarity = correct_hit / len(res_names_1)
        print(self.similarity)
        if self.similarity < 0.4:
            raise ValueError("Sequence Identity less than 40%")

        self.num_residues = len(utilised_res_ind_1)
        print(self.num_residues)
        self.protein_1.utilised_atoms_coords = np.asarray(coords_1)
        self.protein_1.utilised_res_indices = np.asarray(utilised_res_ind_1)
        self.protein_2.utilised_atoms_coords = np.asarray(coords_2)
        self.protein_2.utilised_res_indices = np.asarray(utilised_res_ind_2)
        print("Sequence Checked")

    def check_rmsd(self):
        """
        Superimposes the entire chain of Protein 2 onto Protein 1 using the Ca atoms to see if the RMSD is more than 2.
        """
        ptype = self.protein_1.get_polymer().check_polymer_type()

        chain_superimpose_result: gemmi.SupResult = gemmi.calculate_superposition(self.protein_1.get_polymer(),
                                                                                  self.protein_2.get_polymer(),
                                                                                  ptype, gemmi.SupSelect.CaP)
        print(chain_superimpose_result.rmsd)
        if chain_superimpose_result.rmsd < 2:
            raise ValueError("Proteins RMSD less than 2")

    # def dist_mat_processing(self):
    #     """
    #     Handles the distance matrices of the proteins. First connects to the database to see if the proteins with
    #     the distance matrices exist. If it does not exist, the distance matrix is created and stored into the database.
    #     If something goes wrong when connecting to the database, go straight to offline data management, which is using
    #     the disk storage.
    #     :return:
    #     """
    #     while True:
    #         if conn is None:
    #             self.is_db_connected = False
    #             break
    #         result = check_protein_exists(self.protein_1.code, self.protein_1.chain_param)
    #         if result == -1:
    #             self.is_db_connected = False
    #             break
    #         elif not result:
    #             self.protein_1.get_distance_matrix()
    #             result = insert_protein_dist_mat(self.protein_1.code, self.protein_1.chain_param, self.protein_1.distance_matrix)
    #             if result == -1:
    #                 self.is_db_connected = False
    #                 break
    #         else:
    #             result = get_protein_dist_mat(self.protein_1.code, self.protein_1.chain_param)
    #             if type(result) == int:
    #                 self.is_db_connected = False
    #                 self.protein_1.get_distance_matrix()
    #                 break
    #             self.protein_1.distance_matrix = result
    #
    #         result = check_protein_exists(self.protein_2.code, self.protein_2.chain_param)
    #         if result == -1:
    #             self.is_db_connected = False
    #             break
    #         elif not result:
    #             self.protein_2.get_distance_matrix()
    #             result = insert_protein_dist_mat(self.protein_2.code, self.protein_2.chain_param,
    #                                              self.protein_2.distance_matrix)
    #             if result == -1:
    #                 self.is_db_connected = False
    #                 break
    #         else:
    #             result = get_protein_dist_mat(self.protein_2.code, self.protein_2.chain_param)
    #             if type(result) == int:
    #                 self.is_db_connected = False
    #                 self.protein_2.get_distance_matrix()
    #                 break
    #             self.protein_2.distance_matrix = result
    #         break
    #
    #     if not self.is_db_connected and self.protein_1.distance_matrix is None:
    #         self.protein_1.get_distance_matrix()
    #         self.protein_2.get_distance_matrix()
    #     elif not self.is_db_connected and self.protein_2.distance_matrix is None:
    #         self.protein_2.get_distance_matrix()

    def dist_mat_processing(self):
        """
        Handles the distance matrices of the proteins. First connects to the database to see if the proteins with
        the distance matrices exist. If it does not exist, the distance matrix is created and stored into the database.
        If something goes wrong when connecting to the database, go straight to offline data management, which is using
        the disk storage.
        :return:
        """
        self.protein_1.get_distance_matrix()
        self.protein_2.get_distance_matrix()

    def get_ca_atoms_coords(self):
        """
        Get the coordinates of CA atoms from the proteins. The CA atoms can only be used if the sequence number of the
        residue the atom is in of Protein 1 and 2 are the same. This is to account for protein chains of different lengths.
        The sequence number of residues in the PDB format usually start at 1, but there will be polymer residues in
        protein chains that start at a number above 1 or even a negative value. Returns 4 lists.
        The coordinates of CA atoms and the indices of the residues with utilised CA atoms from proteins 1 and 2.
        The returned lists are the same size.

        Extra residues: CA atoms in residues with sequence numbers 2 and above will be used
        Protein 1 Residue Sequence Numbers [ *  * * * 2 3 ...]
        Protein 2 Residue Sequence Numbers [-2 -1 0 1 2 3 ...]

        Missing residues: CA atom in residues with sequence number 5, 6, 7 will not be used
        Protein 1 Residue Sequence Numbers [ 1 2 3 4 * * * 8 9 ...]
        Protein 2 Residue Sequence Numbers [ 1 2 3 4 5 6 7 8 9 ...]

        :return atom_coords_1: List of coordinates of the CA atoms in Protein 1
        :return atom_coords_2: List of coordinates of the CA atoms in Protein 2
        :return utilised_res_ind_1: List of the indices of the residues that are used in Protein 1. Not the sequence ID number.
        :return utilised_res_ind_2: List of the indices of the residues that are used in Protein 2. Not the sequence ID number.
        """
        # Get the protein 1 and 2 polymer chains. This excludes residues which are only water.
        protein_1_polymer = self.protein_1.get_polymer()
        protein_2_polymer = self.protein_2.get_polymer()
        # Get the length of the polymers
        protein_1_size = len(protein_1_polymer)
        protein_2_size = len(protein_2_polymer)
        # The indices to iterate the polymers
        index_1 = 0
        index_2 = 0
        # Stores the residue names
        res_names_1 = []
        res_names_2 = []
        # Stores the coordinates of the CA atoms
        atom_coords_1 = []
        atom_coords_2 = []
        # Stores the index of the CA atoms in the chains
        utilised_res_ind_1 = []
        utilised_res_ind_2 = []

        if (protein_2_polymer[0].seqid.num - protein_1_polymer[0].seqid.num) > 10:
            seq_num_diff = protein_2_polymer[0].seqid.num - protein_1_polymer[0].seqid.num
        else:
            seq_num_diff = 0

        # Stops iterating at the end of the shorter protein chain
        while index_1 < protein_1_size and index_2 < protein_2_size:
            # Get the sequence id number of the residue in the proteins
            res_num_1 = protein_1_polymer[index_1].seqid.num
            res_num_2 = protein_2_polymer[index_2].seqid.num - seq_num_diff
            # If the residue number of the
            if res_num_1 < res_num_2:
                index_1 += 1
                continue
            if res_num_2 < res_num_1:
                index_2 += 1
                continue
            atom_coord_1 = None
            atom_coord_2 = None
            for a in protein_1_polymer[index_1]:
                if a.name == "CA":
                    atom_coord_1 = a.pos.tolist()
                    break
            for a in protein_2_polymer[index_2]:
                if a.name == "CA":
                    atom_coord_2 = a.pos.tolist()
                    break
            if atom_coord_1 is not None and atom_coord_2 is not None:
                res_names_1.append(protein_1_polymer[index_1].name)
                res_names_2.append(protein_2_polymer[index_2].name)
                utilised_res_ind_1.append(index_1)
                utilised_res_ind_2.append(index_2)
                atom_coords_1.append(atom_coord_1)
                atom_coords_2.append(atom_coord_2)
            index_1 += 1
            index_2 += 1

        return atom_coords_1, atom_coords_2, utilised_res_ind_1, utilised_res_ind_2, res_names_1, res_names_2

    def create_distance_difference_matrix(self):
        """
        Get the distance difference matrix by subtracting one matrix with the other. All values must be positive.
        :return:
        """
        self.diff_dist_mat_init = np.absolute(self.protein_1.distance_matrix - self.protein_2.distance_matrix)
        self.clusters = {i: [i] for i in range(self.diff_dist_mat_init.shape[0])}
        self.link_mat = np.empty((self.diff_dist_mat_init.shape[0] - 1, 4))
        # Save the difference distance matrix image and array before setting the diagonals to infinity for a cleaner visual.
        save_results_to_disk(
            self.output_path,
            self.protein_1.code,
            self.protein_1.chain_param,
            self.protein_2.code,
            self.protein_2.chain_param,
            self.spat_prox,
            self.small_node,
            self.clust_size,
            self.magnitude,
            self.diff_dist_mat_init,
            "diff_dist_mat"
        )
        return self.diff_dist_mat_init

    def run(self):
        # print_diff_dist_mat(self.diff_dist_mat_init)
        np.fill_diagonal(self.diff_dist_mat_init, np.inf)
        start = default_timer()
        diff_dist_mat = np.copy(self.diff_dist_mat_init)
        n = 0
        print("Clustering")
        while len(self.clusters) > 1:
            # print(self.clusters)
            # print(len(self.clusters), self.clusters)
            diff_dist_mat = self.hierarchical_clustering(diff_dist_mat, n)
            if type(diff_dist_mat) == int:
                print("No cluster pair found")
                print(diff_dist_mat)
                print(self.clusters)
                print(len(self.clusters))
                break
            n += 1
        # print("Done")
        end = default_timer()
        total_time = end - start
        # print(self.clusters)
        # print("Time:", total_time)
        save_results_to_disk(
            self.output_path,
            self.protein_1.code,
            self.protein_1.chain_param,
            self.protein_2.code,
            self.protein_2.chain_param,
            self.spat_prox,
            self.small_node,
            self.clust_size,
            self.magnitude,
            self.link_mat,
            "dendrogram"
        )
        superimpose_result = write_to_pdb(self.output_path, self.protein_1, self.protein_2, self.spat_prox, self.small_node, self.clust_size, self.magnitude)
        write_domains_to_pml(self.output_path, self.protein_1, self.protein_2, self.spat_prox, self.small_node, self.clust_size, self.magnitude, self.nodes)
        write_info_file(self.output_path, self.protein_1, self.protein_2, self.spat_prox, self.small_node, self.clust_size, self.magnitude, self.nodes, superimpose_result)
        proteins_str = f"{self.protein_1.code}_{self.protein_1.chain_param}_{self.protein_2.code}_{self.protein_2.chain_param}"
        params_str = f"sp_{self.spat_prox}_node_{self.small_node}_clust_{self.clust_size}_mag_{self.magnitude}"
        print(total_time)
        return round(total_time, 2), len(self.nodes), proteins_str, params_str

    def hierarchical_clustering(self, diff_dist_mat, n):
        """
        Perform hierarchical clustering using the distance difference matrix.
        :param diff_dist_mat:
        :param n: The iteration number
        :return:
        """
        visited_clusters = None
        while True:
            # print("Finding the closest clusters")
            # Find the closest pairs of clusters.
            cluster_pair, min_dist = self.get_closest_clusters(diff_dist_mat, visited_clusters)
            # If no cluster pairs are found, exit early
            if cluster_pair is None:
                print("No cluster pair found")
                return -1
            # print("Checking spatial proximity")
            # When cluster pair is found, check if at least one Ca atom pair between the cluster pair meets the spatial
            # proximity measure.
            spatial_met = self.spatial_proximity_measure(cluster_pair)
            if not spatial_met:
                # print("Spatial proximity not met")
                # If spatial proximity measure is not met, add it to the list of visited clusters. This cluster pair
                # will be ignored in the next iteration when looking for the minimum value.
                if visited_clusters is None:
                    visited_clusters = np.asarray([cluster_pair])
                else:
                    visited_clusters = np.vstack((visited_clusters, cluster_pair))
                continue
            clust_1_size = len(self.clusters[cluster_pair[0]])
            clust_2_size = len(self.clusters[cluster_pair[1]])
            if min_dist >= self.magnitude and clust_1_size > self.small_node and clust_2_size > self.small_node and (clust_1_size + clust_2_size) >= self.clust_size:
                self.clusters[cluster_pair[0]].sort()
                self.clusters[cluster_pair[1]].sort()
                cluster_0_size = len(self.clusters[cluster_pair[0]])
                cluster_1_size = len(self.clusters[cluster_pair[1]])
                # print("Cluster pair distance", cluster_pair, min_dist)
                # print("Cluster pair 0", self.clusters[cluster_pair[0]])
                # print("Cluster pair 1", self.clusters[cluster_pair[1]])
                if cluster_0_size > cluster_1_size:
                    self.nodes[len(self.nodes)] = {
                        "magnitude": min_dist,
                        "large_domain": deepcopy(self.clusters[cluster_pair[0]]),
                        "small_domain": deepcopy(self.clusters[cluster_pair[1]])
                    }
                else:
                    self.nodes[len(self.nodes)] = {
                        "magnitude": min_dist,
                        "large_domain": deepcopy(self.clusters[cluster_pair[1]]),
                        "small_domain": deepcopy(self.clusters[cluster_pair[0]])
                    }
                # self.print_node()
            # print(n, cluster_pair)
            new_cluster_id = n + self.num_residues
            # print(n, self.clusters[cluster_pair[0]], self.clusters[cluster_pair[1]])
            # write_clustering("motion_tree_test", f"{n} {cluster_pair} {min_dist} {self.clusters[cluster_pair[0]]} {self.clusters[cluster_pair[1]]}", n)
            self.clusters[new_cluster_id] = self.clusters[cluster_pair[0]]
            self.clusters[new_cluster_id].extend(self.clusters[cluster_pair[1]])
            # print(cluster_pair[0], cluster_pair[1], min_dist, len(self.clusters[new_cluster_id]))
            self.link_mat[n] = np.array([cluster_pair[0], cluster_pair[1], min_dist, len(self.clusters[new_cluster_id])])
            # print("Update distance matrix")
            # Update the difference distance matrix
            new_diff_dist_mat = self.update_distance_matrix(diff_dist_mat, cluster_pair, new_cluster_id)
            # Remove cluster 1
            del self.clusters[cluster_pair[1]]
            del self.clusters[cluster_pair[0]]
            return new_diff_dist_mat

    def get_closest_clusters(self, diff_dist_matrix: np.array, visited_clusters):
        """
        Using the difference distance matrix, find the closest clusters. If the cluster pair has already been checked,
        ignore it and find another one.
        :param diff_dist_matrix: A NxN matrix of the distance difference matrix
        :param visited_clusters: An array of visited cluster pairs
        :return:
        """
        # Create a copy of the difference distance matrix for dealing with visited clusters
        diff_dist_copy = np.copy(diff_dist_matrix)
        # Set values at the visited cluster indices as infinity so that the next minimum values will be found
        if visited_clusters is not None:
            # If visited clusters are np.array([[3, 5], [7, 10]])
            # The matrix at [3, 5], [5, 3], [7, 10], and [10, 7] will need to be set as infinity.
            rows = np.concatenate((visited_clusters[:, 0], visited_clusters[:, 1]))
            cols = np.concatenate((visited_clusters[:, 1], visited_clusters[:, 0]))
            # https://stackoverflow.com/questions/7761393/how-to-modify-a-2d-numpy-array-at-specific-locations-without-a-loop
            # diff_dist_copy[[3, 7, 5, 10], [5, 10, 3, 7]]
            diff_dist_copy[rows, cols] = np.inf
        # Get the index of the cluster pair with the smallest value.
        # np.argmin gets the index of the minimum value. axis=None means that the array is flattened.
        # np.unravel_index gets the 2D index by specifying the shape of the array
        # https://stackoverflow.com/questions/48135736/what-is-an-intuitive-explanation-of-np-unravel-index
        cluster_pair = np.unravel_index(np.argmin(diff_dist_copy, axis=None), diff_dist_copy.shape)

        return list(cluster_pair), np.min(diff_dist_copy)

    # def spatial_proximity_measure(self, cluster_pair):
    #     """
    #     From the 2 most similar clusters, checks if the closest Ca atom pair meets the spatial proximity
    #     :return:
    #     """
    #     # print(cluster_pair)
    #     cluster_1_indices_list = self.clusters[cluster_pair[0]]
    #     cluster_2_indices_list = self.clusters[cluster_pair[1]]
    #
    #     protein_1_cluster_1_coords = self.protein_1.utilised_atoms_coords[cluster_1_indices_list]
    #     protein_1_cluster_2_coords = self.protein_1.utilised_atoms_coords[cluster_2_indices_list]
    #     protein_2_cluster_1_coords = self.protein_2.utilised_atoms_coords[cluster_1_indices_list]
    #     protein_2_cluster_2_coords = self.protein_2.utilised_atoms_coords[cluster_2_indices_list]
    #
    #     protein_1_distances = cdist(protein_1_cluster_1_coords, protein_1_cluster_2_coords)
    #     protein_2_distances = cdist(protein_2_cluster_1_coords, protein_2_cluster_2_coords)
    #
    #     return np.min(protein_1_distances) < self.spat_prox and np.min(protein_2_distances) < self.spat_prox

    def spatial_proximity_measure(self, cluster_pair):
        """
        From the 2 most similar clusters, checks if the closest Ca atom pair meets the spatial proximity
        :return:
        """
        # print(cluster_pair)
        cluster_1_indices_list = self.clusters[cluster_pair[0]]
        cluster_2_indices_list = self.clusters[cluster_pair[1]]

        is_near_1 = False
        is_near_2 = False

        for i in cluster_1_indices_list:
            for j in cluster_2_indices_list:
                if self.protein_1.distance_matrix[i, j] < self.spat_prox:
                    is_near_1 = True
                if self.protein_2.distance_matrix[i, j] < self.spat_prox:
                    is_near_2 = True

            if is_near_1 and is_near_2:
                break

        return is_near_1 and is_near_2

    def update_distance_matrix(self, diff_dist_mat, cluster_pair, new_id):
        """
        Updates the difference distance matrix after the clustering.
        :param diff_dist_mat: The distance difference matrix
        :param cluster_pair: The pair of cluster IDs which are most similar
        :param new_id: The new cluster id
        :return:
        """
        # https://stackoverflow.com/questions/8486294/how-do-i-add-an-extra-column-to-a-numpy-array
        new_distance_matrix = np.c_[np.copy(diff_dist_mat), np.full(diff_dist_mat.shape[0], np.inf)]
        new_distance_matrix = np.r_[new_distance_matrix, np.full((1, new_distance_matrix.shape[1]), np.inf)]
        new_distance_matrix[cluster_pair, :] = np.inf
        new_distance_matrix[:, cluster_pair] = np.inf
        diss_count = 20
        for k in self.clusters.keys():
            # Ignore the cluster pairs
            if k != cluster_pair[0] and k != cluster_pair[1]:
                # print(self.clusters.keys())
                # print(cluster_pair)
                # print(k, new_id)
                dists = self.get_diff_distances(self.clusters[k], self.clusters[new_id])
                # dists = []
                # dists.extend(self.get_diff_disterences(self.clusters[k], self.clusters[cluster_pair[0]]))
                # dists.extend(self.get_diff_disterences(self.clusters[k], self.clusters[cluster_pair[1]]))
                if len(dists) > diss_count:
                    dists.sort(reverse=True)
                    new_distance_matrix[new_id, k] = mean(dists[:diss_count])
                else:
                    new_distance_matrix[new_id, k] = mean(dists)
                new_distance_matrix[k, new_id] = new_distance_matrix[new_id, k]

        return new_distance_matrix

    def get_diff_distances(self, cluster_indices_1, cluster_indices_2):
        """
        Get the distance differences
        :param cluster_indices_1:
        :param cluster_indices_2:
        :return:
        """
        diff_dists = []
        for i in cluster_indices_1:
            for j in cluster_indices_2:
                diff_dists.append(self.diff_dist_mat_init[i, j])

        return diff_dists

    def print_node(self):
        for key, values in self.nodes.items():
            print(key, "{")
            print("large_domain:", values["large_domain"])
            print("small_domain:", values["small_domain"])
            print("}")


def print_diff_dist_mat(dist_mat):
    row = [str(dist_mat[i]).ljust(4) for i in range(dist_mat.shape[1])]
    print("[".ljust(3), " ", " ".join(row))
    for i in range(dist_mat.shape[0]):
        row = [str(round(j, 1)).ljust(4) for j in dist_mat[i]]
        print(str(i).ljust(3), " ", " ".join(row))
    print("]")


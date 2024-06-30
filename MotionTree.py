from statistics import mean
import numpy as np
from Protein import Protein
from scipy.spatial.distance import cdist
import timeit
from FileMngr import save_results


class MotionTreeTest:
    def __init__(self, files, params=None):
        # Default parameters to be used if not given input
        self.files = files
        if params is None:
            self.parameters = {
                "atoms": "ca",
                "threshold": 7,
                "linkage": "average",
                "avg_diff": 20
            }
        else:
            self.parameters = params
        # Initialise proteins
        self.protein_1: Protein = Protein(
            name=self.files['protein1'],
            file_path=f"{self.files['input_path']}/{self.files['protein1']}.pdb",
            chain=self.files["chain1id"],
            atom_type=self.parameters["atoms"]
        )
        self.protein_2: Protein = Protein(
            name=self.files['protein2'],
            file_path=f"{self.files['input_path']}/{self.files['protein2']}.pdb",
            chain=self.files["chain2id"],
            atom_type=self.parameters["atoms"]
        )
        # The starting distance difference matrix
        self.dist_diff_mat_init = None
        self.num_residues = self.protein_1.main_atoms_coords.shape[0]
        self.create_difference_distance_matrix()
        # Dictionary storing the clusters containing the index of the atoms
        self.clusters = {i: [i] for i in range(self.dist_diff_mat_init.shape[0])}
        # The linkage matrix for building the dendrogram
        self.link_mat = np.empty((self.dist_diff_mat_init.shape[0] - 1, 4))
        # Bending residues indices
        self.bending_residues = set()

    def run(self):
        start = timeit.default_timer()
        dist_diff_mat = np.copy(self.dist_diff_mat_init)
        n = 0
        while len(self.clusters) > 1:
            # print(self.clusters)
            dist_diff_mat = self.hierarchical_clustering(dist_diff_mat, n)
            # print(dist_diff_mat.shape)
            # print(dist_diff_mat)
            # print(self.clusters)
            if type(dist_diff_mat) == int:
                print("No cluster pair found")
                print(dist_diff_mat)
                print(self.clusters)
                print(len(self.clusters))
                break
            n += 1

        end = timeit.default_timer()
        # print(self.clusters)
        print("Time:", end - start)
        save_results(
            self.files["output_path"],
            f"{self.protein_1.name}_{self.protein_1.chain_param}_{self.protein_2.name}_{self.protein_2.chain_param}",
            self.link_mat,
            "dendrogram"
        )

    def create_difference_distance_matrix(self):
        """
        Get the distance difference matrix by subtracting one matrix with the other. All values must be positive.
        :return:
        """
        self.dist_diff_mat_init = np.absolute(self.protein_1.distance_matrix - self.protein_2.distance_matrix)
        # Save the difference distance matrix before setting the diagonals to infinity for a cleaner visual.
        save_results(
            self.files["output_path"],
            f"{self.protein_1.name}_{self.protein_1.chain_param}_{self.protein_2.name}_{self.protein_2.chain_param}",
            self.dist_diff_mat_init,
            "diff_dist_mat"
        )
        np.fill_diagonal(self.dist_diff_mat_init, np.inf)

    def hierarchical_clustering(self, dist_diff_mat, n):
        """
        Perform hierarchical clustering using the distance difference matrix.
        :param dist_diff_mat:
        :param n: The iteration number
        :return:
        """
        visited_clusters = None
        while True:
            # print("Finding the closest clusters")
            # Find the closest pairs of clusters.
            cluster_pair, min_dist = self.get_closest_clusters(dist_diff_mat, visited_clusters)
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
            new_dist_diff_mat = self.update_distance_matrix(dist_diff_mat, cluster_pair, new_cluster_id)
            # Remove cluster 1
            del self.clusters[cluster_pair[1]]
            del self.clusters[cluster_pair[0]]
            return new_dist_diff_mat

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

        # print(cluster_pair, np.min(diff_dist_copy))

        return list(cluster_pair), np.min(diff_dist_copy)

    def spatial_proximity_measure(self, cluster_pair):
        """
        From the 2 most similar clusters, checks if the closest Ca atom pair meets the spatial proximity
        :return:
        """
        cluster_1_indices_list = self.clusters[cluster_pair[0]]
        cluster_2_indices_list = self.clusters[cluster_pair[1]]

        protein_1_cluster_1_coords = self.protein_1.main_atoms_coords[cluster_1_indices_list]
        protein_1_cluster_2_coords = self.protein_1.main_atoms_coords[cluster_2_indices_list]
        protein_2_cluster_1_coords = self.protein_2.main_atoms_coords[cluster_1_indices_list]
        protein_2_cluster_2_coords = self.protein_2.main_atoms_coords[cluster_2_indices_list]

        protein_1_distances = cdist(protein_1_cluster_1_coords, protein_1_cluster_2_coords)
        protein_2_distances = cdist(protein_2_cluster_1_coords, protein_2_cluster_2_coords)

        return np.min(protein_1_distances) < self.parameters["threshold"] and np.min(protein_2_distances) < self.parameters["threshold"]

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
        # new_distance_matrix[:, cluster_pair[0]] = np.inf
        # new_distance_matrix[cluster_pair[0], :] = np.inf
        # new_distance_matrix[:, cluster_pair[1]] = np.inf
        # new_distance_matrix[cluster_pair[1], :] = np.inf
        new_distance_matrix[cluster_pair, :] = np.inf
        new_distance_matrix[:, cluster_pair] = np.inf
        # if cluster_pair[0] == 287:
        #     print(new_distance_matrix[cluster_pair[0], :])
        for k in self.clusters.keys():
            # Ignore the cluster pairs
            if k != cluster_pair[0] and k != cluster_pair[1]:
                # print(self.clusters.keys())
                # print(cluster_pair)
                # print(k, new_id)
                dists = self.get_dist_differences(self.clusters[k], self.clusters[new_id])
                # dists = []
                # dists.extend(self.get_dist_differences(self.clusters[k], self.clusters[cluster_pair[0]]))
                # dists.extend(self.get_dist_differences(self.clusters[k], self.clusters[cluster_pair[1]]))
                if len(dists) > 20:
                    dists.sort(reverse=True)
                    new_distance_matrix[new_id, k] = mean(dists[:20])
                else:
                    new_distance_matrix[new_id, k] = mean(dists)
                new_distance_matrix[k, new_id] = new_distance_matrix[new_id, k]

        return new_distance_matrix

    def get_dist_differences(self, cluster_indices_1, cluster_indices_2):
        """
        Get the distance differences
        :param cluster_indices_1:
        :param cluster_indices_2:
        :return:
        """
        dist_diffs = []
        for i in cluster_indices_1:
            for j in cluster_indices_2:
                dist_diffs.append(self.dist_diff_mat_init[i, j])

        return dist_diffs


def print_dist_diff_mat(dist_mat):
    for i in range(dist_mat.shape[0]):
        row = []
        for j in range(dist_mat.shape[1]):
            row.append(dist_mat[i, j])
        print(row)


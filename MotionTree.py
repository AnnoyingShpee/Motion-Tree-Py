import sys
from statistics import mean
import numpy as np
import pandas as pd
from Protein import Protein
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform, cdist
import matplotlib.pyplot as plt
import timeit

class MotionTree:
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
            file_path=f"{self.files['input_path']}{self.files['protein1']}.pdb",
            chain=self.files["chain1id"],
            atom_type=self.parameters["atoms"]
        )
        self.protein_2: Protein = Protein(
            name=self.files['protein2'],
            file_path=f"{self.files['input_path']}{self.files['protein2']}.pdb",
            chain=self.files["chain2id"],
            atom_type=self.parameters["atoms"]
        )
        # The starting distance difference matrix
        self.dist_diff_mat_init = None
        self.create_difference_distance_matrix()
        # Dictionary storing the clusters containing the index of the atoms
        self.clusters = {i: [i] for i in range(self.dist_diff_mat_init.shape[0])}
        # The linkage matrix for building the dendrogram
        self.linkage_mat = np.array([])
        # Bending residues indices
        self.bending_residues = set()

    def run(self):
        # dist_diff_mat = self.first_hierarchical_clustering()
        # print(self.clusters)
        # if type(dist_diff_mat) == int:
        #     print("No cluster pair found.")
        #     sys.exit()
        link_mat = np.empty((self.dist_diff_mat_init.shape[0] - 1, 4))
        forest = {i: [i] for i in range(self.dist_diff_mat_init.shape[0])}
        start = timeit.default_timer()
        dist_diff_mat = np.copy(self.dist_diff_mat_init)
        n = 0
        while len(self.clusters) > 1:
            dist_diff_mat, link, forest = self.hierarchical_clustering(dist_diff_mat, forest)
            link_mat[n] = link
            n += 1
            # print(self.clusters)
            if type(dist_diff_mat) == int:
                print("No cluster pair found")
                print(dist_diff_mat)
                print(self.clusters)
                print(len(self.clusters))
                break

        end = timeit.default_timer()
        # print(self.clusters)
        print("Time:", end - start)
        # self.plot_clusters(self.protein_1.main_atoms_coords, self.clusters)

        # print("Initial")
        # self.print_dist_diff_mat(self.dist_diff_mat_init)

        # print("Final")
        # self.print_dist_diff_mat(dist_diff_mat)

        # condensed_dist_diff = squareform(1-dist_diff_mat, force="tovector", checks=False)
        # print(condensed_dist_diff)
        # print(condensed_dist_diff.shape)

        # sch.dendrogram(sch.linkage(squareform(dist_diff_mat, force="tovector", checks=False)))
        # sch.dendrogram(linkage_mat)
        # plt.show()
        plt.matshow(self.dist_diff_mat_init)
        # plt.savefig("Init.png")
        # plt.matshow(dist_diff_mat)
        # plt.show()
        # plt.savefig("Final.png")
        # plt.savefig(f"data/output/dendrograms/{self.protein_1.name}_{self.protein_1.chain_param}_{self.protein_2.name}_{self.protein_2.chain_param}_dendro.png")
        # plt.matshow(dist_diff_mat)
        # plt.savefig(f"data/output/matrices/{self.protein_1.name}_{self.protein_1.chain_param}_{self.protein_2.name}_{self.protein_2.chain_param}_mat.png")

    def create_difference_distance_matrix(self):
        """
        Get the distance difference matrix by subtracting one matrix with the other. All values must be positive.
        :return:
        """
        self.dist_diff_mat_init = np.absolute(self.protein_1.distance_matrix - self.protein_2.distance_matrix)
        # Save the difference distance matrix before setting the diagonals to infinity for a cleaner visual.
        plt.matshow(self.dist_diff_mat_init)
        plt.savefig(f"data/output/matrices/{self.protein_1.name}_{self.protein_1.chain_param}_{self.protein_2.name}_{self.protein_2.chain_param}_mat.png")
        np.fill_diagonal(self.dist_diff_mat_init, np.inf)

    def hierarchical_clustering(self, dist_diff_mat, forest):
        """
        Perform hierarchical clustering using the distance difference matrix.
        :param dist_diff_mat:
        :param forest:
        :return:
        """
        visited_clusters = None
        while True:
            # print("Finding the closest clusters")
            # Find the closest pairs of clusters.
            cluster_pair = self.get_closest_clusters(dist_diff_mat, visited_clusters)
            # If no cluster pairs are found, exit early
            if cluster_pair is None:
                print("No cluster pair found")
                return -1
            # print("Checking spatial proximity")
            # When cluster pair is found, check if at least one Ca atom pair between the cluster pair meets the spatial
            # proximity measure.
            if not self.spatial_proximity_measure(cluster_pair):
                # print("Spatial proximity not met")
                # If spatial proximity measure is not met, add it to the list of visited clusters. This cluster pair
                # will be ignored in the next iteration when looking for the minimum value
                if visited_clusters is None:
                    visited_clusters = np.asarray([cluster_pair])
                else:
                    visited_clusters = np.vstack((visited_clusters, cluster_pair))
                continue

            self.clusters[cluster_pair[0]].extend(self.clusters[cluster_pair[1]])
            # print("Update distance matrix")
            new_dist_diff_mat = self.update_distance_matrix(dist_diff_mat, cluster_pair)
            del self.clusters[cluster_pair[1]]
            return new_dist_diff_mat

    def get_closest_clusters(self, diff_dist_matrix: np.array, visited_clusters):
        """
        Using the difference distance matrix, find the closest clusters. If the cluster pair has already been checked,
        ignore it and find another one.
        :param diff_dist_matrix: A NxN matrix of the distance difference matrix
        :param visited_clusters: An array of visited cluster pairs
        :return:
        """
        # Create a copy of the matrix for dealing with visited clusters
        diff_dist_copy = np.copy(diff_dist_matrix)
        # Set values at the visited cluster indices as infinity so that the next minimum values will be found
        if visited_clusters is not None:
            rows = np.concatenate((visited_clusters[:, 0], visited_clusters[:, 1]))
            cols = np.concatenate((visited_clusters[:, 1], visited_clusters[:, 0]))
            diff_dist_copy[rows, cols] = np.inf
        # print(np.unravel_index(np.argmin(diff_dist_copy, axis=None), diff_dist_copy.shape))
        cluster_pair = np.unravel_index(np.argmin(diff_dist_copy, axis=None), diff_dist_copy.shape)

        return list(cluster_pair)

    # def get_closest_clusters(self, diff_dist_matrix: np.array, visited_clusters):
    #     """
    #     Using the difference distance matrix, find the closest clusters. If the cluster pair has already been checked,
    #     ignore it and find another one.
    #     :param diff_dist_matrix:
    #     :param visited_clusters:
    #     :return:
    #     """
    #     diff_dist_copy = np.copy(diff_dist_matrix)
    #     np.fill_diagonal(diff_dist_copy, np.inf)
    #
    #     visited_clusters_set = set(visited_clusters)
    #     upper_diag_indices = np.triu_indices(diff_dist_matrix.shape[0], k=1)
    #     print(upper_diag_indices)
    #     print(*upper_diag_indices)
    #     print(zip(*upper_diag_indices))
    #     # print(visited_clusters)
    #     if visited_clusters is not None:
    #         # print(visited_clusters)
    #         # print(diff_dist_copy[visited_clusters])
    #         rows = np.concatenate((visited_clusters[:, 0], visited_clusters[:, 1]))
    #         cols = np.concatenate((visited_clusters[:, 1], visited_clusters[:, 0]))
    #         diff_dist_copy[rows, cols] = np.inf
    #     print(np.unravel_index(np.argmin(diff_dist_copy, axis=None), diff_dist_copy.shape))
    #     cluster_pair = np.unravel_index(np.argmin(diff_dist_copy, axis=None), diff_dist_copy.shape)
    #
    #     return list(cluster_pair)

    def spatial_proximity_measure(self, cluster_pair):
        """
        From the 2 most similar clusters, checks if the Ca atoms closest to each other meet the spatial proximity
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

        if protein_1_distances.shape[0] > 1:
            np.fill_diagonal(protein_1_distances, np.inf)
            np.fill_diagonal(protein_2_distances, np.inf)

        return np.min(protein_1_distances) < self.parameters["threshold"] and np.min(protein_2_distances) < self.parameters["threshold"]

    def update_distance_matrix(self, diff_dist_mat, cluster_pair):
        """

        :param diff_dist_mat: The distance difference matrix
        :param cluster_pair: The pair of cluster IDs which are most similar
        :return:
        """
        new_distance_matrix = np.copy(diff_dist_mat)
        for k in self.clusters.keys():
            # Ignore the cluster pairs
            if k != cluster_pair[0] and k != cluster_pair[1]:
                dists = []
                dists.extend(self.get_avg_dist_differences(self.clusters[k], self.clusters[cluster_pair[0]]))
                dists.extend(self.get_avg_dist_differences(self.clusters[k], self.clusters[cluster_pair[1]]))
                if len(dists) > 20:
                    dists.sort(reverse=True)
                    new_distance_matrix[cluster_pair[0], k] = mean(dists[:20])
                else:
                    new_distance_matrix[cluster_pair[0], k] = mean(dists)
                new_distance_matrix[k, cluster_pair[0]] = new_distance_matrix[cluster_pair[0], k]
        new_distance_matrix[:, cluster_pair[1]] = np.inf
        new_distance_matrix[cluster_pair[1], :] = np.inf
        return new_distance_matrix

    def get_avg_dist_differences(self, cluster_indices_1, cluster_indices_2):
        dist_diffs = []

        for i in cluster_indices_1:
            for j in cluster_indices_2:
                dist_diffs.append(self.dist_diff_mat_init[i, j])

        return dist_diffs

    def print_dist_diff_mat(self, dist_mat):
        for i in range(dist_mat.shape[0]):
            row = []
            for j in range(dist_mat.shape[1]):
                row.append(dist_mat[i, j])
            print(row)


# def solution(x, y):
#     a = list(zip(x, y))  # This produces list of tuples
#     ax = sorted(a, key=lambda x: x[0])  # Presorting x-wise
#     ay = sorted(a, key=lambda x: x[1])  # Presorting y-wise
#     p1, p2, mi = closest_pair(ax, ay)  # Recursive D&C function
#     return mi
#
#
# def closest_pair(ax, ay):
#     ln_ax = len(ax)  # It's quicker to assign variable
#     if ln_ax <= 3:
#         return brute(ax)  # A call to bruteforce comparison
#     mid = ln_ax // 2  # Division without remainder, need int
#     Qx = ax[:mid]  # Two-part split
#     Rx = ax[mid:]
#     # Determine midpoint on x-axis
#     midpoint = ax[mid][0]
#     Qy = list()
#     Ry = list()
#     for x in ay:  # split ay into 2 arrays using midpoint
#         if x[0] <= midpoint:
#            Qy.append(x)
#         else:
#            Ry.append(x)
#     # Call recursively both arrays after split
#     (p1, q1, mi1) = closest_pair(Qx, Qy)
#     (p2, q2, mi2) = closest_pair(Rx, Ry)
#     # Determine smaller distance between points of 2 arrays
#     if mi1 <= mi2:
#         d = mi1
#         mn = (p1, q1)
#     else:
#         d = mi2
#         mn = (p2, q2)
#     # Call function to account for points on the boundary
#     (p3, q3, mi3) = closest_split_pair(ax, ay, d, mn)
#     # Determine smallest distance for the array
#     if d <= mi3:
#         return mn[0], mn[1], d
#     else:
#         return p3, q3, mi3
#
#
# def brute(ax):
#     mi = dist(ax[0], ax[1])
#     p1 = ax[0]
#     p2 = ax[1]
#     ln_ax = len(ax)
#     if ln_ax == 2:
#         return p1, p2, mi
#     for i in range(ln_ax-1):
#         for j in range(i + 1, ln_ax):
#             if i != 0 and j != 1:
#                 d = dist(ax[i], ax[j])
#                 if d < mi:  # Update min_dist and points
#                     mi = d
#                     p1, p2 = ax[i], ax[j]
#     return p1, p2, mi
#
# def closest_split_pair(p_x, p_y, delta, best_pair):
#     ln_x = len(p_x)  # store length - quicker
#     mx_x = p_x[ln_x // 2][0]  # select midpoint on x-sorted array
#     # Create a subarray of points not further than delta from
#     # midpoint on x-sorted array
#     s_y = [x for x in p_y if mx_x - delta <= x[0] <= mx_x + delta]
#     best = delta  # assign best value to delta
#     ln_y = len(s_y)  # store length of subarray for quickness
#     for i in range(ln_y - 1):
#         for j in range(i+1, min(i + 7, ln_y)):
#             p, q = s_y[i], s_y[j]
#             dst = dist(p, q)
#             if dst < best:
#                 best_pair = p, q
#                 best = dst
#     return best_pair[0], best_pair[1], best
#

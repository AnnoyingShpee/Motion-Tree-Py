import sys
from statistics import mean
import numpy as np
import pandas as pd
from Protein import Protein
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
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
        self.dist_diff_mat_init = self.create_difference_distance_matrix()
        # Dictionary storing the clusters containing the index of the atoms
        self.clusters = {i: [i] for i in range(self.dist_diff_mat_init.shape[0])}
        # Bending residues indices
        self.bending_residues = set()

    def run(self):
        # dist_diff_mat = self.first_hierarchical_clustering()
        # print(self.clusters)
        # if type(dist_diff_mat) == int:
        #     print("No cluster pair found.")
        #     sys.exit()

        dist_diff_mat = np.copy(self.dist_diff_mat_init)
        while len(self.clusters) > 1:
            start = timeit.default_timer()
            dist_diff_mat = self.hierarchical_clustering(dist_diff_mat)
            end = timeit.default_timer()
            print("Time:", end - start)
            print(self.clusters)
            if type(dist_diff_mat) == int:
                print("No cluster pair found")
                print(dist_diff_mat)
                print(self.clusters)
                print(len(self.clusters))
                break

        # self.plot_clusters(self.protein_1.main_atoms_coords, self.clusters)

        # print("Initial")
        # self.print_dist_diff_mat(self.dist_diff_mat_init)

        # print("Final")
        # self.print_dist_diff_mat(dist_diff_mat)

        # condensed_dist_diff = squareform(1-dist_diff_mat, force="tovector", checks=False)
        # print(condensed_dist_diff)
        # print(condensed_dist_diff.shape)
        sch.dendrogram(sch.linkage(squareform(dist_diff_mat, force="tovector", checks=False)))
        plt.savefig(f"data/output/dendrograms/{self.protein_1.name}_{self.protein_1.chain_param}_{self.protein_2.name}_{self.protein_2.chain_param}_dendro.png")
        plt.matshow(dist_diff_mat)
        plt.savefig(f"data/output/matrices/{self.protein_1.name}_{self.protein_1.chain_param}_{self.protein_2.name}_{self.protein_2.chain_param}_mat.png")

    def create_difference_distance_matrix(self):
        return np.absolute(self.protein_1.distance_matrix - self.protein_2.distance_matrix)

    # def first_hierarchical_clustering(self):
    #     new_dist_diff_mat = None
    #     visited_clusters = []
    #     while new_dist_diff_mat is None:
    #         cluster_pair = self.closest_clusters(self.dist_diff_mat_init, visited_clusters)
    #         if cluster_pair is None:
    #             return -1
    #         visited_clusters.append(cluster_pair)
    #         if not self.spatial_proximity_measure(cluster_pair):
    #             print("Spatial proximity not met")
    #             continue
    #         self.clusters[cluster_pair[0]].extend(self.clusters[cluster_pair[1]])
    #         new_dist_diff_mat = self.update_distance_matrix(self.dist_diff_mat_init, cluster_pair)
    #         del self.clusters[cluster_pair[1]]
    #
    #     return new_dist_diff_mat

    def hierarchical_clustering(self, dist_diff_mat):
        visited_clusters = []
        cluster_found = False
        while not cluster_found:
            # print("Finding closest clusters")
            cluster_pair = self.closest_clusters(dist_diff_mat, visited_clusters)
            if cluster_pair is None:
                print("No cluster pair found")
                return -1
            visited_clusters.append(cluster_pair)
            # print("Checking spatial proximity")
            if not self.spatial_proximity_measure(cluster_pair):
                # print("Spatial proximity not met")
                continue
            self.clusters[cluster_pair[0]].extend(self.clusters[cluster_pair[1]])
            # print("Update distance matrix")
            new_dist_diff_mat = self.update_distance_matrix(dist_diff_mat, cluster_pair)
            del self.clusters[cluster_pair[1]]
            return new_dist_diff_mat

    def closest_clusters(self, diff_dist_matrix: np.array, visited_clusters):
        min_distance = np.inf
        cluster_pair = None
        for i in self.clusters:
            for j in self.clusters:
                if i != j and diff_dist_matrix[i, j] < min_distance and (i, j) not in visited_clusters and (j, i) not in visited_clusters:
                    min_distance = diff_dist_matrix[i, j]
                    cluster_pair = (i, j)

        return cluster_pair

    def spatial_proximity_measure(self, cluster_pair):
        """
        From the 2 most similar clusters, checks if the Ca atoms closest to each other meet the spatial proximity
        :return:
        """
        cluster_1_indices_list = self.clusters[cluster_pair[0]]
        cluster_2_indices_list = self.clusters[cluster_pair[1]]
        min_distance_1 = np.inf
        min_distance_2 = np.inf
        for i in cluster_1_indices_list:
            for j in cluster_2_indices_list:
                if self.protein_1.distance_matrix[i, j] < min_distance_1:
                    min_distance_1 = self.protein_1.distance_matrix[i, j]
                if self.protein_2.distance_matrix[i, j] < min_distance_2:
                    min_distance_2 = self.protein_2.distance_matrix[i, j]

        return min_distance_1 < self.parameters["threshold"] and min_distance_2 < self.parameters["threshold"]

    def update_distance_matrix(self, diff_dist_mat, cluster_pair):
        """

        :param diff_dist_mat: The distance difference matrix
        :param cluster_pair: The pair of cluster IDs which are most similar
        :return:
        """
        new_distance_matrix = np.copy(diff_dist_mat)
        for k in self.clusters.keys():
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
        new_distance_matrix[cluster_pair[0], cluster_pair[0]] = 0
        return new_distance_matrix

    def get_avg_dist_differences(self, cluster_indices_1, cluster_indices_2):
        dist_diffs = []

        for i in cluster_indices_1:
            for j in cluster_indices_2:
                dist_diffs.append(self.dist_diff_mat_init[i, j])

        return dist_diffs

    def plot_clusters(self, chain, clusters):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        for cluster, points in clusters.items():
            for point in points:
                ax.scatter(chain[point][0], chain[point][1], chain[point][2])
        plt.show()

    def print_dist_diff_mat(self, dist_mat):
        for i in range(dist_mat.shape[0]):
            row = []
            for j in range(dist_mat.shape[1]):
                row.append(dist_mat[i, j])
            print(row)


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


import sys
from statistics import mean
import numpy as np
import pandas as pd
from Protein import Protein
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform, cdist
import matplotlib.pyplot as plt
import timeit


class MotionTreeInit:
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
        start = timeit.default_timer()
        dist_diff_mat = np.copy(self.dist_diff_mat_init)
        n = 0
        while len(self.clusters) > 1:
            dist_diff_mat = self.hierarchical_clustering(dist_diff_mat, n)
            # print(self.clusters)
            if type(dist_diff_mat) == int:
                print("No cluster pair found")
                print(dist_diff_mat)
                print(self.clusters)
                print(len(self.clusters))
                break
            n += 1
        end = timeit.default_timer()
        print(self.clusters)
        # print(np.unique(self.clusters[0], return_counts=True))
        print("Time:", end - start)
        return self.clusters

    def create_difference_distance_matrix(self):
        """
        Get the distance difference matrix by subtracting one matrix with the other. All values must be positive.
        :return:
        """
        self.dist_diff_mat_init = np.absolute(self.protein_1.distance_matrix - self.protein_2.distance_matrix)
        # Save the difference distance matrix before setting the diagonals to infinity for a cleaner visual.
        plt.matshow(self.dist_diff_mat_init)
        plt.savefig(
            f"data/output/matrices/{self.protein_1.name}_{self.protein_1.chain_param}_{self.protein_2.name}_{self.protein_2.chain_param}_mat.png")
        np.fill_diagonal(self.dist_diff_mat_init, np.inf)

    def hierarchical_clustering(self, dist_diff_mat, n):
        """
        Perform hierarchical clustering using the distance difference matrix.
        :param dist_diff_mat:
        :return:
        """
        visited_clusters = []
        cluster_found = False
        while not cluster_found:
            # print("Finding the closest clusters")
            cluster_pair, min_dist = self.get_closest_clusters(dist_diff_mat, visited_clusters, n)
            if cluster_pair is None:
                print("No cluster pair found")
                return -1
            # if n > 230:
            #     print(cluster_pair)
            # print("Checking spatial proximity")
            spatial_met, clust_1, clust_2, min_dist_1, min_dist_2 = self.spatial_proximity_measure(cluster_pair)
            if not spatial_met:
                # print("Spatial proximity not met")
                visited_clusters.append(cluster_pair)
                continue
            # print(n, cluster_pair)
            self.clusters[cluster_pair[0]].extend(self.clusters[cluster_pair[1]])
            # print("Update distance matrix")
            new_dist_diff_mat = self.update_distance_matrix(dist_diff_mat, cluster_pair)
            del self.clusters[cluster_pair[1]]
            return new_dist_diff_mat

    def get_closest_clusters(self, diff_dist_matrix: np.array, visited_clusters, n):
        """
        Using the difference distance matrix, find the closest clusters. If the cluster pair has already been checked,
        ignore it and find another one.
        :param diff_dist_matrix:
        :param visited_clusters:
        :return:
        """
        min_distance = np.inf
        cluster_pair = None
        for i in self.clusters:
            for j in self.clusters:
                if i != j and diff_dist_matrix[i, j] < min_distance and (i, j) not in visited_clusters and (j, i) not in visited_clusters:
                    min_distance = diff_dist_matrix[i, j]
                    cluster_pair = (i, j)

        return cluster_pair, min_distance

    def spatial_proximity_measure(self, cluster_pair):
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


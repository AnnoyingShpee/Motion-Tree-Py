import sys
from statistics import mean
import numpy as np
import pandas as pd
from Protein import Protein
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt


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
        dist_diff_mat = self.first_hierarchical_clustering()
        # print(self.clusters)
        if type(dist_diff_mat) == int:
            print("No cluster pair found.")
            sys.exit()

        while len(self.clusters) > 1:
            dist_diff_mat = self.hierarchical_clustering(dist_diff_mat)
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
        sch.dendrogram(sch.linkage(dist_diff_mat))
        # plt.matshow(dist_diff_mat)
        plt.show()

    def create_difference_distance_matrix(self):
        return np.absolute(self.protein_1.distance_matrix - self.protein_2.distance_matrix)

    def first_hierarchical_clustering(self):
        new_dist_diff_mat = None
        visited_clusters = []
        while new_dist_diff_mat is None:
            cluster_pair = self.closest_clusters(self.dist_diff_mat_init, visited_clusters)
            if cluster_pair is None:
                return -1
            visited_clusters.append(cluster_pair)
            if not self.spatial_proximity_measure(cluster_pair):
                print("Spatial proximity not met")
                continue
            self.clusters[cluster_pair[0]].extend(self.clusters[cluster_pair[1]])
            new_dist_diff_mat = self.update_distance_matrix(self.dist_diff_mat_init, cluster_pair)
            del self.clusters[cluster_pair[1]]

        return new_dist_diff_mat

    def hierarchical_clustering(self, dist_diff_mat):
        new_dist_diff_mat = dist_diff_mat
        visited_clusters = []
        cluster_found = False
        while not cluster_found:
            cluster_pair = self.closest_clusters(new_dist_diff_mat, visited_clusters)
            if cluster_pair is None:
                print("No cluster pair found")
                return -1
            visited_clusters.append(cluster_pair)
            if not self.spatial_proximity_measure(cluster_pair):
                # print("Spatial proximity not met")
                continue
            self.clusters[cluster_pair[0]].extend(self.clusters[cluster_pair[1]])
            cluster_found = True
            new_dist_diff_mat = self.update_distance_matrix(new_dist_diff_mat, cluster_pair)
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

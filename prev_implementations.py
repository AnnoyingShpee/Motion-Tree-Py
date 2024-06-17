# def hierarchical_clustering(self, dist_diff_mat, forest):
#     """
#     Perform hierarchical clustering using the distance difference matrix.
#     :param dist_diff_mat:
#     :param forest:
#     :return:
#     """
#     visited_clusters = []
#     # visited_clusters = None
#     cluster_found = False
#     while not cluster_found:
#         # print("Finding the closest clusters")
#         cluster_pair = self.get_closest_clusters(dist_diff_mat, visited_clusters)
#         if cluster_pair is None:
#             print("No cluster pair found")
#             return -1
#         visited_clusters.append(cluster_pair)
#         # print(cluster_pair)
#         # if visited_clusters is None:
#         #     visited_clusters = np.asarray([cluster_pair])
#         # else:
#         #     visited_clusters = np.vstack((visited_clusters, cluster_pair))
#         # print("Checking spatial proximity")
#         if not self.spatial_proximity_measure(cluster_pair):
#             # print("Spatial proximity not met")
#             continue
#         self.clusters[cluster_pair[0]].extend(self.clusters[cluster_pair[1]])
#         # print("Update distance matrix")
#         new_dist_diff_mat = self.update_distance_matrix(dist_diff_mat, cluster_pair)
#         del self.clusters[cluster_pair[1]]
#         return new_dist_diff_mat


# def get_closest_clusters(self, diff_dist_matrix: np.array, visited_clusters):
#     """
#     Using the difference distance matrix, find the closest clusters. If the cluster pair has already been checked,
#     ignore it and find another one.
#     :param diff_dist_matrix:
#     :param visited_clusters:
#     :return:
#     """
#     min_distance = np.inf
#     cluster_pair = None
#     for i in self.clusters:
#         for j in self.clusters:
#             if i != j and diff_dist_matrix[i, j] < min_distance and (i, j) not in visited_clusters and (j, i) not in visited_clusters:
#                 # print(i, j, diff_dist_matrix[i, j])
#                 min_distance = diff_dist_matrix[i, j]
#                 cluster_pair = (i, j)
#
#     return cluster_pair


# def spatial_proximity_measure(self, cluster_pair):
#     cluster_1_indices_list = self.clusters[cluster_pair[0]]
#     cluster_2_indices_list = self.clusters[cluster_pair[1]]
#     min_distance_1 = np.inf
#     min_distance_2 = np.inf
#     for i in cluster_1_indices_list:
#         for j in cluster_2_indices_list:
#             if self.protein_1.distance_matrix[i, j] < min_distance_1:
#                 min_distance_1 = self.protein_1.distance_matrix[i, j]
#             if self.protein_2.distance_matrix[i, j] < min_distance_2:
#                 min_distance_2 = self.protein_2.distance_matrix[i, j]
#
#     return min_distance_1 < self.parameters["threshold"] and min_distance_2 < self.parameters["threshold"]


# def update_distance_matrix(self, diff_dist_mat, cluster_pair):
#     """
#
#     :param diff_dist_mat: The distance difference matrix
#     :param cluster_pair: The pair of cluster IDs which are most similar
#     :return:
#     """
#     new_distance_matrix = np.copy(diff_dist_mat)
#     for k in self.clusters.keys():
#         if k != cluster_pair[0] and k != cluster_pair[1]:
#             dists = []
#             dists.extend(self.get_avg_dist_differences(self.clusters[k], self.clusters[cluster_pair[0]]))
#             dists.extend(self.get_avg_dist_differences(self.clusters[k], self.clusters[cluster_pair[1]]))
#             if len(dists) > 20:
#                 dists.sort(reverse=True)
#                 new_distance_matrix[cluster_pair[0], k] = mean(dists[:20])
#             else:
#                 new_distance_matrix[cluster_pair[0], k] = mean(dists)
#             new_distance_matrix[k, cluster_pair[0]] = new_distance_matrix[cluster_pair[0], k]
#     return new_distance_matrix
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
from symnmf import initialize_decomposition_matrix_H
import symnmfmodule as c
import math
import sys
np.random.seed(1234)

EPSILON = 0.0001


class Cluster:
    """
    A class to represent a cluster from HW1
    """
    def __init__(self, centroid: list):
        """
        Initialize a new cluster
        :param centroid: the centroid of the cluster
        """
        self.sum: list = [0.0 for i in range(len(centroid))]
        """
        sum of vector's coordinates
        """
        self.centroid: list = centroid
        self.size = 0
        """
        temp size for each rotation
        """

    def get_size(self):
        """
        Get the size of the cluster
        """
        return self.size

    def get_sum(self):
        """
        Get the sum of values in the cluster
        """
        return self.sum

    def get_centroid(self):
        """
        Get the centroid of the cluster
        """
        return self.centroid

    def set_size(self, new_size):
        self.size = new_size

    def set_sum(self, new_sum):
        self.sum = new_sum

    def add_xi(self, vect_xi):
        """
        Add a vector to the cluster
        :param vect_xi: The vector to add.
        """
        for i in range(len(vect_xi)):
            self.sum[i] += vect_xi[i]
        self.size += 1

    def update_centroid(self):
        """
        Update to the new centroid of the cluster
        """
        new_centroid = [0 for i in range(len(self.centroid))]
        for i in range(len(self.sum)):
            new_centroid[i] = float(self.sum[i]) / self.size
        self.centroid = new_centroid

    def set_centroid(self, new_centroid):
        self.centroid = new_centroid

    def reset_sum_and_size(self):
        """
        Reset cluster's sum and size for a new loop.
        """
        self.sum = [0 for i in range(len(self.centroid))]
        self.size = 0

    def calculate_euclidean_distance(self, vect_xi):
        """
        Calculate euclidean distance between a new vector and the cluster's centroid
        :param vect_xi: The vector we wish to calculate the euclidean distance from the cluster's centroid.
        :return:
        """
        sum = 0
        for i in range(len(vect_xi)):
            sum += math.pow(vect_xi[i] - self.centroid[i], 2)
        return math.sqrt(sum)

    def __repr__(self):
        return str(self.centroid)


def initialize_centroids(vect_arr, K):
    """
    Initialize all clusters with K initial vectors from data points as centroids.
    :param vect_arr: The data points
    :param K: Number of clusters.
    :return: The initialized clusters
    """
    clusters = [Cluster(vect_arr[i]) for i in range(K)]
    return clusters


def calculate_closest_cluster(clusters, vect_xi):
    """
    Calculate the closest cluster to the given vector.
    :param clusters: The list of clusters.
    :param vect_xi: The vector we're working with.
    :return: The cluster closest to the given vector and its index.
    """
    min_eucledian_distance = sys.maxsize
    min_cluster = None
    min_cluster_index = 0
    for i, cluster in enumerate(clusters):
        temp_ed = cluster.calculate_euclidean_distance(vect_xi)
        if temp_ed < min_eucledian_distance:
            min_eucledian_distance = temp_ed
            min_cluster = cluster
            min_cluster_index = i
    return min_cluster, min_cluster_index


def k_means_algorithm(vectors, K, iter_limit=300):
    """
    Runs the K-Means algorithm from HW1.
    :param vectors: Data points.
    :param K: Number of clusters.
    :param iter_limit: The maximum number of iterations
    :return: A list of the assignment of all vectors to their cluster and a list of final centroids.
    """
    clusters = initialize_centroids(vectors, K)
    """ Initialize starting K clusters """
    labels = [-1] * len(vectors)
    iter_number = 0
    flag = False
    while iter_number <= iter_limit and not flag:
        iter_number += 1
        for i in range(len(vectors)):
            xi = vectors[i]
            """ Assign every xi to the closest cluster k """
            min_cluster, min_cluster_index = calculate_closest_cluster(clusters, xi)
            min_cluster.add_xi(xi)
            labels[i] = min_cluster_index
        flag = True
        """ Update centroids """
        for cluster in clusters:
            """ Update the centroids and check for convergence """
            prev_cluster_centroid = cluster.get_centroid()
            cluster.update_centroid()
            if flag:
                convergence = cluster.calculate_euclidean_distance(prev_cluster_centroid)
                if convergence >= EPSILON:
                    flag = False
            cluster.reset_sum_and_size()
    for cluster in clusters:
        cluster.set_centroid(["%.4f" % xi for xi in cluster.get_centroid()])
    centroids = [cluster.get_centroid() for cluster in clusters]
    return labels, centroids


def load_data(file_name: str) -> list:
    """
    Load dataset from a .txt file
    :param file_name: The name of the .txt file.
    :return: A numpy array representing the data.
    """
    data_points = pd.read_csv(file_name, header=None)
    data_points = data_points.to_numpy()
    return data_points.tolist()


def apply_kmeans(data_points, k: int):
    """
    Apply KMeans clustering
    :param data_points: Data vectors represented in a 2D array.
    :param k: number of clusters.
    :return: KMeans clustering of all vectors.
    """
    labels, centroids = k_means_algorithm(data_points, k, 300)
    return labels


def symnmf_clustering(H):
    """
    Derive clusters by finding the index of the max value in each row of H
    :param H: The decomposition matrix calculated in step 1.4
    :return: A clustering solution for the vectors
    """
    return np.argmax(H, axis=1)


def apply_symnmf(data_points, k: int):
    """
    Apply SymNMF clustering
    :param data_points: Data vectors represented in a 2D array.
    :param k: number of clusters.
    :return: SymNMF clustering of all vectors.
    """
    vectors_count, vector_length = len(data_points), len(data_points[0])
    normalized_similarity_matrix = c.norm(data_points)
    if normalized_similarity_matrix is None:
        return
    m = np.mean(normalized_similarity_matrix)
    initial_H = initialize_decomposition_matrix_H(vectors_count, m, k)
    result_matrix = c.symnmf(k, initial_H, normalized_similarity_matrix, vectors_count)
    if result_matrix is None:
        return
    labels = symnmf_clustering(result_matrix)
    return labels


def get_silhouette_score(data_points, labels):
    """
    Calculate Silhouette Score
    :param data_points: Data vectors represented in a 2D array
    :param labels:
    :return: The Silhouette Score
    """
    return silhouette_score(data_points, labels)


def main():
    """
    Main function
    """
    file_name = sys.argv[2]
    k = int(sys.argv[1])
    data_points = load_data(file_name)

    # SymNMF clustering
    symnmf_labels = apply_symnmf(data_points, k)
    symnmf_silhouette = get_silhouette_score(data_points, symnmf_labels)
    print(f"nmf: {symnmf_silhouette:.4f}")

    # KMeans clustering
    kmeans_labels = apply_kmeans(data_points, k)
    kmeans_silhouette = get_silhouette_score(data_points, kmeans_labels)
    print(f"kmeans: {kmeans_silhouette:.4f}")


if __name__ == "__main__":
    main()

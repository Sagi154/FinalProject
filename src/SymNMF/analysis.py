import sys
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import symnmf
from kmeans import *
from src.SymNMF.symnmf import parse_arguments
import symnmfmodule as c
np.random.seed(1234)
import math

# Load dataset from .txt file
def load_data(file_name):
    return np.loadtxt(file_name)



def avg_W_entries(normalized_similarity_matrix, vectors_count):
    total_matrix_sum = 0
    # total_matrix_sum = np.sum(np.array(normalized_similarity_matrix))
    for line in normalized_similarity_matrix:
        for element in line:
            total_matrix_sum += element
    return total_matrix_sum / (vectors_count * vectors_count)



def initialize_decomposition_matrix_H(vectors_count, m, K):
    # Upper bound not tight
    decomposition_matrix = [[np.random.uniform(low=0, high=(2 * math.sqrt(m/K) + 1e-10)) for j in range(K)] for i in range(vectors_count)]
    return decomposition_matrix



# Apply KMeans clustering
def apply_kmeans(data_points, k):
    labels, centroids = k_means_algorithm(data_points, k, 200)
    return labels


def symnmf_clustering(H):
    # Derive clusters by finding the index of the max value in each row of H
    return np.argmax(H, axis=1)



# Apply SymNMF clustering
def apply_symnmf(data_points, k):
    vectors_count, vector_length = len(data_points), len(data_points[0])
    normalized_similarity_matrix = c.norm(data_points, vectors_count, vector_length)
    if normalized_similarity_matrix is None:
        return
    m = avg_W_entries(normalized_similarity_matrix, vectors_count)
    initial_H = initialize_decomposition_matrix_H(vectors_count, m, k)
    result_matrix = c.symnmf(k, initial_H, normalized_similarity_matrix, vectors_count)
    if result_matrix is None:
        return
    labels = symnmf_clustering(result_matrix)
    return labels


# Calculate Silhouette Score
def get_silhouette_score(data_points, labels):
    return silhouette_score(data_points, labels)


# Main function
def main(file_name, k):
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
    import sys

    file_name = sys.argv[2]
    k = int(sys.argv[1])
    main(file_name, k)

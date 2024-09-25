import sys

from src.SymNMF.symnmf import parse_arguments

def parse_arguments():
    K = int(float(sys.argv[1]))
    file_name = sys.argv[2]
    return K, file_name

def symnmf_cluster_assignment():


def main():
    K, file_name = parse_arguments()


import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import symnmf  # Import your C module for SymNMF
from kmeans import *


# Load dataset from .txt file
def load_data(file_name):
    return np.loadtxt(file_name)


# Apply KMeans clustering
def apply_kmeans(data, k):
    kmeans = KMeans(n_clusters=k, init='random', n_init=10, tol=0.001, max_iter=300)
    labels = kmeans.fit_predict(data)
    return labels


def symnmf_clustering(H):
    # Derive clusters by finding the index of the max value in each row of H
    return np.argmax(H, axis=1)



# Apply SymNMF clustering
def apply_symnmf(data, k):
    # Generate the similarity matrix A, degree matrix D, and normalized similarity matrix W
    A = symnmf.sym(data)  # Assuming this computes the similarity matrix A
    D = symnmf.ddg(A)  # Assuming this computes the degree matrix D
    W = symnmf.norm(A, D)  # Assuming this normalizes similarity to get W

    # Factorize W to get H using SymNMF
    H = symnmf.symnmf(W, k)  # SymNMF factorization returns matrix H

    # Assign clusters based on the maximum value in each row of H
    labels = symnmf_clustering(H)
    return labels


# Calculate Silhouette Score
def get_silhouette_score(data, labels):
    return silhouette_score(data, labels)


# Main function
def main(file_name, k):
    data = load_data(file_name)

    # KMeans clustering
    kmeans_labels = apply_kmeans(data, k)
    kmeans_silhouette = get_silhouette_score(data, kmeans_labels)
    print(f"KMeans Silhouette Score: {kmeans_silhouette:.4f}")

    # SymNMF clustering
    symnmf_labels = apply_symnmf(data, k)
    symnmf_silhouette = get_silhouette_score(data, symnmf_labels)
    print(f"SymNMF Silhouette Score: {symnmf_silhouette:.4f}")


if __name__ == "__main__":
    import sys

    file_name = sys.argv[2]
    k = int(sys.argv[1])
    main(file_name, k)

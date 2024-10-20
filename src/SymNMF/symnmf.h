# include <stdio.h>
# include <stdlib.h>
# include <math.h>

/** The number of vectors in our data */
extern int vectors_count;

/** The number of coordinates in every vector in our data */
extern int vector_length;

/** The number of clusters */
extern int K;

/** The data vectors */
extern double **data_vectors;

/**
 * @brief Frees the memory of a given matrix.
 *
 * This function frees each (successfully created) row of the matrix individually,
 * and then frees the matrix itself.
 * @param matrix A pointer to the matrix (double**) whose memory is to be freed.
 * @param number_of_rows The number of rows to free, corresponding to the number of successfully created rows.
 * @return None
 */
void free_memory_of_matrix(double **matrix, int number_of_rows);

/**
 * @brief Calculates the similarity matrix.
 *
 * This function calculates and returns the similarity matrix as described in step 1.1.
 * @param None
 * @return A pointer to the similarity matrix (double**).
 */
double **calculate_similarity_matrix();


/**
 * @brief Calculates the diagonal degree matrix.
 *
 * This function calculates and returns the diagonal degree matrix based on the
 * provided similarity matrix as described in step 1.2.
 * @param similarity_matrix A pointer to the similarity matrix (double**) calculated in step 1.1.
 * @return A pointer to the diagonal degree matrix (double**).
 */
double **calculate_diagonal_degree_matrix(double **similarity_matrix);

/**
 * @brief Calculates the normalized similarity matrix.
 *
 * This function calculates and returns the normalized similarity matrix based on the provided
 * diagonal degree matrix and similarity matrix, as described in step 1.3.
 * @param diagonal_degree_matrix A pointer to the diagonal degree matrix (double**) calculated in step 1.2.
 * @param sym_matrix A pointer to the similarity matrix (double**) calculated in step 1.1.
 * @return A pointer to the normalized similarity matrix (double**).
 */
double **calculate_normalized_similarity_matrix(double **diagonal_degree_matrix , double **sym_matrix);


/**
 * @brief Performs the full symNMF algorithm and outputs H.
 *
 * This function performs the symNMF algorithm, updating and returning the decomposition matrix H until convergence.
 * @param decomposition_matrix_H A pointer to the initial decomposition matrix (double**) calculated in Python.
 * @param normalized_similarity_matrix A pointer to the normalized similarity matrix (double**) calculated in step 1.3.
 * @return A pointer to the final decomposition matrix H (double**), once convergence is reached.
 */
double **calculate_final_decomposition_matrix_symnmf(double **decomposition_matrix_H, double **normalized_similarity_matrix);


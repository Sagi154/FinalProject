# include <stdio.h>
# include <stdlib.h>
# include <math.h>

extern int vectors_count;
extern int vector_length;
extern int K;
extern double **data_vectors;


void free_memory_of_matrix(double **matrix, int numb_of_rows);

/*
 * Calculate and output the similarity matrix as described in 1.1.
 */
double **calculate_similarity_matrix();

/*
 * Calculate and output the Diagonal Degree Matrix as described in 1.2.
 */
double **calculate_diagonal_degree_matrix(double ** similarity_matrix);

/*
 * Calculate and output the Diagonal Degree Matrix as described in 1.2.
 */
double **calculate_normalized_similarity_matrix(double** diagonal_degree_matrix , double** sym_matrix);

/*
 * Perform full the symNMF as described in 1 and output H.
 */
double **calculate_final_decomposition_matrix_symnmf(double **decomposition_matrix_H, double **normalized_similarity_matrix);


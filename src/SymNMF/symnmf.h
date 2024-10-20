# include <stdio.h>
# include <stdlib.h>
# include <math.h>

/* The number of vectors in our data */
extern int vectors_count;

/* The number of coordinates in every vector in our data */
extern int vector_length;

/* The number of clusters */
extern int K;

/* The data vectors */
extern double **data_vectors;

/*
Frees the memory of a given matrix by freeing each row's array one by one starting from
the last row (of successfully created rows) and then free the array itself.

Input:
    double **matrix: The matrix whose memory we want to free
    int numb_of_rows: The number of rows we wish to free which is equal to the number of rows
    successfully created.

Output: None
*/
void free_memory_of_matrix(double **matrix, int numb_of_rows);

/*
Calculates and outputs the similarity matrix as described step in 1.1.

Input: None

Output:
    double **sym_matrix: returns the similarity matrix.
*/
double **calculate_similarity_matrix();

/*
Calculates and outputs the Diagonal Degree Matrix as described in step 1.2.

Input:
    double **similarity_matrix: the similarity matrix calculated in step 1.1.

Output:
    double **diagonal_degree_matrix: returns the diagonal degree matrix.
*/
double **calculate_diagonal_degree_matrix(double ** similarity_matrix);

/*
Calculates and outputs the normalized similarity matrix as described in step 1.3.

Input:
    double **diagonal_degree_matrix: the diagonal degree matrix calculated in step 1.2.
    double **sym_matrix: the similarity matrix calculated in step 1.1.

Output:
    double **graph_Laplacian: returns the normalized similarity matrix.
*/
double **calculate_normalized_similarity_matrix(double **diagonal_degree_matrix , double **sym_matrix);


/*
Performs the full symNMF algorithm and outputs H.

Input:
    double **decomposition_matrix_H: the initial decomposition matrix calculated in python.
    double **normalized_similarity_matrix: the normalized similarity matrix calculated in step 1.3.

Output:
    double **decomposition_matrix_H: returns the final decomposition_matrix_H once reached convergence.
*/
double **calculate_final_decomposition_matrix_symnmf(double **decomposition_matrix_H, double **normalized_similarity_matrix);


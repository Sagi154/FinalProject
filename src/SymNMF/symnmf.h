# include <stdio.h>
# include <stdlib.h>
# include <math.h>

struct cluster
{
 int vect_count;
 int vect_length;
 double* centroid;
 double* sum;
} ;

struct cord
{
 double value;
 struct cord *next;
};

struct vector
{
 struct vector *next;
 struct cord *cords;
};

double calculate_squared_euclidean_distance(double *first_vector, double *second_vector);

double calculate_frobenius_norm_squared(double **matrix, int rows_count, int columns_count);

void copy_vector_by_cord(double *copy_from, double *copy_to);

void free_vector_cords(struct cord *head_cord, int cords_counted);

void free_memory_of_lists(struct vector *head_vec, int vectors_counted);

void free_memory_of_matrix(double **matrix, int index);

void free_memory_of_vectors_array(int vectors_counted);

int create_vector_arr(struct vector *head_vec);

int is_double_integer(double value);

double **calculate_inverse_square_root(double **diagonal_matrix);

double **multiply_matrices(double **first_matrix, double **second_matrix, int rows_count, int columns_count, int common_count);

double **transpose_matrix(double **matrix, int rows_count, int columns_count);

double **matrices_subtraction(double** first_matrix, double** second_matrix, int rows_count, int columns_count);

double avg_W_entries(double ** normalized_similarity_matrix);

void initialize_H(double **decomposition_matrix, double m, int rows_count, int columns_count);

double **sym();
/*
Calculate and output the similarity matrix as described in 1.1.
 */
double **ddg(double ** similarity_matrix);
/*
Calculate and output the Diagonal Degree Matrix as described in 1.2
 */
double **norm(double** diagonal_degree_matrix , double** sym_matrix);
/*
Calculate and output the normalized similarity matrix as described in 1.3
 */

double **symnmf(double **normalized_similarity_matrix);
 /*
 Perform full the symNMF as described in 1 and output H.
 */

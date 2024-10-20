#include <assert.h>
#include "symnmf.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define ERR_MSG "An Error Has Occurred"
#define DBL_MAX 1.7976931348623157e+308
#define BETA 0.5
#define EPSILON 0.0001
#define MAX_ITER 300
#define MAX_LINE_LENGTH 1024

int vectors_count, vector_length, K;
double **data_vectors;


/**
 * @brief Calculates the squared Euclidean distance between two vectors.
 *
 * This function computes the squared Euclidean distance between the given two vectors.
 * @param first_vector A pointer to the first vector (double*).
 * @param second_vector A pointer to the second vector (double*).
 * @return The squared Euclidean distance as a double.
 */
double calculate_squared_euclidean_distance(double *first_vector, double *second_vector)
{
    int i;
    double sum_of_points = 0.0;
    for (i = 0 ; i < vector_length; i++)
    {
        sum_of_points += pow(first_vector[i] - second_vector[i], 2);
    }
    return sum_of_points;
}


/**
 * @brief Calculates the squared Frobenius norm of a matrix.
 *
 * This function computes the squared Frobenius norm of the given matrix.
 * @param matrix A pointer to the matrix (double**) for which to calculate the norm.
 * @param rows_count The number of rows in the matrix.
 * @param columns_count The number of columns in the matrix.
 * @return The squared Frobenius norm as a double.
 */
double calculate_frobenius_norm_squared(double **matrix, int rows_count, int columns_count)
{
    int i,j;
    double sum = 0.0;
    for (i = 0; i < rows_count; i++) {
        for (j = 0; j < columns_count; j++) {
            sum += pow(fabs(matrix[i][j]), 2);
        }
    }
    return sum;
}


void free_memory_of_matrix(double **matrix, int number_of_rows)
{
    int i;
    for(i = 0; i < number_of_rows; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}


/**
 * @brief Calculates the inverse square root of a diagonal matrix.
 *
 * This function computes the inverse square root of the provided diagonal matrix.
 * It creates a new matrix such that for every entry on the matrix's diagonal the value
 * is (1/ sqrt(diagonal matrix received, in the same entry)) and every other entry is 0.
 * @param diagonal_matrix A pointer to the diagonal matrix (double**) for which to calculate the inverse square root.
 * @return A pointer to the resulting matrix (double**), which is dynamically allocated.
 */
double **calculate_inverse_square_root(double **diagonal_matrix) {
    int i;
    double **inverse_square_root_matrix = (double **) calloc(vectors_count, sizeof(double*));
    if (inverse_square_root_matrix == NULL) {
        return NULL;
    }
    for (i = 0; i < vectors_count; i++) {
        inverse_square_root_matrix[i] = (double *)calloc(vectors_count, sizeof(double));
        if (inverse_square_root_matrix[i] == NULL)
        {
            free_memory_of_matrix(inverse_square_root_matrix, i);
            return NULL;
        }
        inverse_square_root_matrix[i][i] = 1/ sqrt(diagonal_matrix[i][i]);
    }
    return inverse_square_root_matrix;
}


/**
 * @brief Multiplies two matrices.
 *
 * This function computes the product of the two input matrices.
 * @param first_matrix A pointer to the first matrix (double**) to multiply.
 * @param second_matrix A pointer to the second matrix (double**) to multiply.
 * @param rows_count The number of rows in the resulting matrix.
 * @param columns_count The number of columns in the resulting matrix.
 * @param common_count The number of columns in the first matrix (and rows in the second matrix).
 * @return A pointer to the resulting matrix (double**), which is dynamically allocated.
 */
double **multiply_matrices(double **first_matrix, double **second_matrix, int rows_count, int columns_count, int common_count) {
    int i,j,k;
    double **product_matrix = (double **) calloc(rows_count, sizeof(double *));
    if (product_matrix == NULL) {
        return NULL;
    }
    for (i = 0; i < rows_count; i++) {
        product_matrix[i] = (double *)calloc(columns_count, sizeof(double));
        if (product_matrix[i] == NULL)
        {
            free_memory_of_matrix(product_matrix, i);
            return NULL;
        }
        for (j = 0; j < columns_count; j++) {
            for (k = 0; k < common_count; k++) {
                product_matrix[i][j] += first_matrix[i][k] * second_matrix[k][j];
            }
        }
    }
    return product_matrix;
}


/**
 * @brief Transposes a given matrix.
 *
 * This function computes the transpose of the provided matrix.
 * @param matrix A pointer to the matrix (double**) to be transposed.
 * @param rows_count The number of rows in the resulting matrix.
 * @param columns_count The number of columns in the resulting matrix.
 * @return A pointer to the transposed matrix (double**), which is dynamically allocated.
 */
double **transpose_matrix(double **matrix, int rows_count, int columns_count) {
    int i,j;
    double **transpose_matrix = (double **) calloc(rows_count, sizeof(double *));
    if (transpose_matrix == NULL)
    {
        return NULL;
    }
    for (i = 0; i < rows_count; i++) {
        transpose_matrix[i] = (double *)calloc(columns_count, sizeof(double));
        if (transpose_matrix[i] == NULL)
        {
            free_memory_of_matrix(transpose_matrix, i);
            return NULL;
        }
        for (j = 0; j < columns_count; j++) {
            transpose_matrix[i][j] = matrix[j][i];
        }
    }
    return transpose_matrix;
}


/**
 * @brief Computes the subtraction of two matrices.
 *
 * This function subtracts the second matrix from the first matrix element-wise
 * and returns the resulting matrix.
 * @param first_matrix A pointer to the first matrix (double**) from which to subtract.
 * @param second_matrix A pointer to the second matrix (double**) to be subtracted.
 * @param rows_count The number of rows in the matrices.
 * @param columns_count The number of columns in the matrices.
 * @return A pointer to the resulting matrix (double**), which is dynamically allocated.
 */
double **matrices_subtraction(double** first_matrix, double** second_matrix, int rows_count, int columns_count) {
    int i,j;
    double **subtraction_matrix = (double **) calloc(rows_count, sizeof(double *));
    if (subtraction_matrix == NULL)
    {
        return NULL;
    }
    for (i = 0; i < rows_count; i++) {
        subtraction_matrix[i] = (double *)calloc(columns_count, sizeof(double));
        if (subtraction_matrix[i] == NULL)
        {
            free_memory_of_matrix(subtraction_matrix, i);
            return NULL;
        }
        for (j = 0; j < columns_count; j++) {
            subtraction_matrix[i][j] = first_matrix[i][j] - second_matrix[i][j];
        }
    }
    return subtraction_matrix;
}


/**
 * @brief Prints a matrix to the standard output.
 *
 * This function outputs the elements of the given matrix in the required format.
 * @param matrix A pointer to the matrix (double**) to be printed.
 * @param rows_count The number of rows in the matrix.
 * @param columns_count The number of columns in the matrix.
 * @return None
 */
void print_matrix(double** matrix, int rows_count, int columns_count)
{
    int i, j;
    for (i = 0; i < rows_count; i++)
    {
        for (j = 0; j < columns_count - 1; j++)
        {
            printf("%.4f,", matrix[i][j]);
        }
        printf("%.4f\n", matrix[i][columns_count-1]);
    }
}


double **calculate_similarity_matrix(){
    int i,j;
    double **sym_matrix;
    sym_matrix = (double**) calloc(vectors_count, sizeof(double*));
    if (sym_matrix == NULL)
        return NULL;
    for (i = 0; i < vectors_count; i++){
        sym_matrix[i] = (double *)calloc( vectors_count, sizeof(double));
        if (sym_matrix[i] == NULL)
        {
            free_memory_of_matrix(sym_matrix, i);
            return NULL;
        }
        for (j = 0; j < vectors_count; j++){
            /*
             * Calculates and inserts the value of the similarity matrix at entry (i, j)
             * using the "calculate_squared_euclidean_distance" function according to the formula.
             */
            if(i == j)
                sym_matrix[i][j] = 0.0;
            else
                sym_matrix[i][j] = exp(-0.5 * calculate_squared_euclidean_distance(data_vectors[i], data_vectors[j]));
        }
    }
    return sym_matrix;
}


double **calculate_diagonal_degree_matrix(double **similarity_matrix) {
    int i,j;
    double i_row_sum = 0.0;
    double **diagonal_degree_matrix;
    diagonal_degree_matrix = (double**) calloc(vectors_count, sizeof(double*));
    if (diagonal_degree_matrix == NULL)
        return NULL;
    for(i = 0; i < vectors_count; i++) {
        diagonal_degree_matrix[i] = (double *)calloc(vectors_count, sizeof(double));
        i_row_sum = 0.0;
        if (diagonal_degree_matrix[i] == NULL) {
            free_memory_of_matrix(diagonal_degree_matrix, i);
            return NULL;
        }
        /* Calculates the sum of the row */
        for(j = 0; j < vectors_count; j++) {
            i_row_sum += similarity_matrix[i][j];
        }
        diagonal_degree_matrix[i][i] = i_row_sum;
    }
    return diagonal_degree_matrix;
}


double **calculate_normalized_similarity_matrix(double **diagonal_degree_matrix , double **sym_matrix){
    double **inverse_square_root_matrix, **graph_Laplacian_left, **graph_Laplacian;
    /* Calculates the inverse square root matrix of the diagonal degree matrix that received */
    inverse_square_root_matrix = calculate_inverse_square_root(diagonal_degree_matrix);
    if (inverse_square_root_matrix == NULL)
        return NULL;
    /*
     * Multiplies the inverse square root matrix of the diagonal degree matrix provided with the sym matrix provided
     * such that the inverse square root matrix of the diagonal degree matrix is in the left and the sym matrix in the right
    */
    graph_Laplacian_left = multiply_matrices(inverse_square_root_matrix, sym_matrix , vectors_count, vectors_count, vectors_count);
    if (graph_Laplacian_left == NULL) {
        free_memory_of_matrix(inverse_square_root_matrix, vectors_count);
        return NULL;
    }
    /*
     * Multiplies the product matrix calculated above with the inverse square root matrix of the diagonal degree matrix provided
     * such that the product matrix calculated above is in the left and the inverse square root matrix of the diagonal degree matrix in the right
     */
    graph_Laplacian = multiply_matrices(graph_Laplacian_left, inverse_square_root_matrix, vectors_count, vectors_count, vectors_count);
    if (graph_Laplacian == NULL) {
        free_memory_of_matrix(inverse_square_root_matrix, vectors_count);
        free_memory_of_matrix(graph_Laplacian_left, vectors_count);
        return NULL;
    }
    free_memory_of_matrix(graph_Laplacian_left, vectors_count);
    free_memory_of_matrix(inverse_square_root_matrix, vectors_count);
    return graph_Laplacian;
}


/**
 * @brief Calculates the matrix resulting from the product of 'H', 'H^t', and 'H'.
 *
 * This function computes the matrix 'H*H^t*H' where 'H' is the decomposition matrix.
 * It first transposes the matrix 'H', then computes the product 'H*H^t', and finally
 * multiplies that result by 'H'.
 * @param decomposition_matrix_H A pointer to the decomposition matrix H.
 * @return A pointer to the resulting matrix 'H*H^t*H', or NULL if an error occurs during
 *         matrix operations (e.g., memory allocation failure).
 */
double **calculate_H_Ht_H_matrix(double **decomposition_matrix_H ) {
    double **matrix_H_transposed, **matrix_H_H_t, **matrix_H_Ht_H;
    matrix_H_transposed = transpose_matrix(decomposition_matrix_H, K, vectors_count);
    if (matrix_H_transposed == NULL)
        return NULL;
    matrix_H_H_t = multiply_matrices(decomposition_matrix_H, matrix_H_transposed, vectors_count, vectors_count, K);
    if (matrix_H_H_t == NULL) {
        free_memory_of_matrix(matrix_H_transposed, K);
        return NULL;
    }
    matrix_H_Ht_H = multiply_matrices(matrix_H_H_t, decomposition_matrix_H, vectors_count, K, vectors_count);
    free_memory_of_matrix(matrix_H_transposed, K);
    free_memory_of_matrix(matrix_H_H_t, vectors_count);
    return matrix_H_Ht_H;
}


/**
 * @brief Calculates the updated matrix H based on the given decomposition matrix 'H',
 *        the product matrix 'WH', and the matrix 'H*H^t*H'.
 *
 * This function computes the updated values for the decomposition matrix H using the
 * formula given in the project instructions at step 1.4.2.
 * @param decomposition_matrix_H A pointer to the current decomposition matrix 'H'.
 * @param matrix_WH A pointer to the product matrix 'WH'.
 * @param matrix_H_Ht_H A pointer to the matrix resulting from 'H*H^t*H'.
 * @return A pointer to the updated matrix H, or NULL if an error occurs during memory allocation
 *         or matrix operations.
 */
double **calculate_updated_H(double **decomposition_matrix_H, double **matrix_WH, double **matrix_H_Ht_H) {
    int i, j;
    double **updated_H;
    updated_H = (double**)calloc(vectors_count, sizeof(double*));
    if (updated_H == NULL) {
        free_memory_of_matrix(matrix_WH, vectors_count);
        free_memory_of_matrix(matrix_H_Ht_H, vectors_count);
        return NULL;
    }
    for (i = 0; i < vectors_count; i++) {
        updated_H[i] = (double *)calloc(K, sizeof(double));
        if (updated_H[i] == NULL) {
            free_memory_of_matrix(matrix_WH, vectors_count);
            free_memory_of_matrix(matrix_H_Ht_H, vectors_count);
            free_memory_of_matrix(updated_H, i);
            return NULL;
        }
        for (j = 0; j < K; j++) {
            updated_H[i][j] = decomposition_matrix_H[i][j] * (1 - BETA + BETA * (matrix_WH[i][j] / matrix_H_Ht_H[i][j]));
        }
    }
    free_memory_of_matrix(matrix_WH, vectors_count);
    free_memory_of_matrix(matrix_H_Ht_H, vectors_count);
    return updated_H;
}


double **calculate_final_decomposition_matrix_symnmf(double **decomposition_matrix_H, double **normalized_similarity_matrix) {
    int iter_count;
    double frobenius_norm_value;
    double **matrix_WH, **matrix_H_Ht_H, **subtraction_matrix, **updated_H;
    iter_count = 1;
    while (iter_count <= MAX_ITER) {
        matrix_WH = multiply_matrices(normalized_similarity_matrix, decomposition_matrix_H, vectors_count, K, vectors_count);
        if (matrix_WH == NULL) {
            free_memory_of_matrix(decomposition_matrix_H, vectors_count);
            return NULL; }
        matrix_H_Ht_H = calculate_H_Ht_H_matrix(decomposition_matrix_H);
        if (matrix_H_Ht_H == NULL) {
            free_memory_of_matrix(decomposition_matrix_H, vectors_count);
            free_memory_of_matrix(matrix_WH, vectors_count); }
        updated_H = calculate_updated_H(decomposition_matrix_H, matrix_WH, matrix_H_Ht_H);
        if(updated_H == NULL) {
            free_memory_of_matrix(decomposition_matrix_H, vectors_count);
            return NULL; }
        /* Calculates the required matrix for convergence check */
        subtraction_matrix = matrices_subtraction(updated_H, decomposition_matrix_H , vectors_count , K);
        if (subtraction_matrix == NULL) {
            free_memory_of_matrix(decomposition_matrix_H, vectors_count);
            free_memory_of_matrix(updated_H, vectors_count);
            return NULL; }
        frobenius_norm_value = calculate_frobenius_norm_squared(subtraction_matrix, vectors_count, K);
        free_memory_of_matrix(subtraction_matrix, vectors_count);
        free_memory_of_matrix(decomposition_matrix_H, vectors_count);
        decomposition_matrix_H = updated_H;
        /* Checks for convergence */
        if (frobenius_norm_value < EPSILON)
            break;
        iter_count++;
    }
    return decomposition_matrix_H;
}


/**
 * @brief Retrieves the dimensions of a matrix from a file.
 *
 * This function reads the first line of the specified file to determine the number of
 * elements in each vector (vector length) and counts the total number of vectors in the file.
 * It updates the global variables `vectors_count` and `vector_length` accordingly.
 * @param vectors_file A pointer to the file (FILE*) containing the matrix data.
 * @return 0 if the dimensions are successfully retrieved, or 1 if there is an error.
 */
/* TODO: only returns 0? */
int get_matrix_dimensions(FILE* vectors_file) {
    char line[MAX_LINE_LENGTH];
    vectors_count = 0;
    vector_length = 0;
    if (fgets(line, MAX_LINE_LENGTH, vectors_file)) {
        /*
         * Reads the first vector in order to find how many rows there are in each vector and
         * does that by counting these rows one by one
         */
        char *token = strtok(line, ",");
        while (token != NULL) {
        vector_length++;
        token = strtok(NULL, ",");
        }
    }
    vectors_count++;
    while(fgets(line, MAX_LINE_LENGTH, vectors_file)) {
        /* Counts how many vectors there are in the given input starting from the second one (the first one counted above) */
        vectors_count++;
    }

    /* Returns to the beginning of the file in order to read the vectors afterward */
    rewind(vectors_file);
    return 0;
}


/**
 * @brief Reads vector data from a file and stores it in a dynamically allocated array.
 *
 * This function reads vectors from the specified file, allocating memory for each
 * vector and storing the coordinates as double values.
 * @param vectors_file A pointer to the file (FILE*) containing the vector data.
 * @return 0 if the vectors are successfully read, 1 if there is an error.
 */
int read_vectors(FILE* vectors_file) {
    int j, i;
    char line[MAX_LINE_LENGTH];
    char *token;
    i = 0;
    data_vectors = (double**) calloc(vectors_count, sizeof(double*));
    if(data_vectors == NULL)
        return 1;
    while(fgets(line, MAX_LINE_LENGTH, vectors_file)) {
        data_vectors[i] = (double*) calloc(vector_length, sizeof(double));
        if(data_vectors[i] == NULL) {
            free_memory_of_matrix(data_vectors, i);
            return 1;
        }
        /* Reads the values of the coordinates of vector i using the "," delimiter */
        token = strtok(line, ",");
        for(j = 0; j < vector_length; j++) {
            if(token != NULL) {
                data_vectors[i][j] = atof(token);
                token = strtok(NULL, ",");
            }
        }
        i++;
    }
    return 0;
}


/**
 * @brief Extracts matrix data from a specified file.
 *
 * This function opens the provided file, retrieves the matrix dimensions using
 * `get_matrix_dimensions`, and reads the vector data using `read_vectors`.
 * @param file_name A pointer to the file name (FILE*) to read the matrix data from.
 * @return 0 if the data is successfully extracted, or 1 if there is an error.
 */
int extract_data_from_file(FILE* file_name) {
    FILE* vectors_file;
    vectors_file = fopen(file_name, "r");
    if(vectors_file == NULL)
        return 1;
    get_matrix_dimensions(vectors_file);
    if(read_vectors(vectors_file) != 0)
        return 1;
    fclose(vectors_file);
    return 0;
}


/**
 * @brief Performs the specified goal operation on vector data.\n \n
 * This function calculates the similarity matrix, diagonal degree matrix, or normalized similarity
 * matrix based on the provided goal string. It prints the resulting matrix if applicable.
 * @param goal A string indicating the goal operation: "sym" for similarity matrix,
 *             "ddg" for diagonal degree matrix, or "norm" for normalized similarity matrix.
 * @return 0 if the operation is successfully performed, or 1 if there is an error.
 */
int perform_goal(char* goal) {
    double **sym_matrix, **ddg_matrix, **normalized_similarity_matrix;
    /* Calculates the similarity matrix of the vectors array matrix, the method will use it for any goal. */
    sym_matrix = calculate_similarity_matrix();
    if(sym_matrix == NULL)
        return 1;
    if(strcmp(goal, "sym") == 0)
        /* If the goal is to calculate the similarity matrix: print it and finish. */
        print_matrix(sym_matrix, vectors_count, vectors_count);
    else {
        /* If the goal isn't to calculate the similarity matrix, the method calculates the diagonal degree matrix
         * of  similarity matrix calculated above, the method will now use it for any goal. */
        ddg_matrix = calculate_diagonal_degree_matrix(sym_matrix);
        if(ddg_matrix == NULL) {
            free_memory_of_matrix(sym_matrix, vectors_count);
            return 1; }
        if(strcmp(goal, "ddg") == 0)
            /* If the goal is to calculate the diagonal degree matrix: print it and finish. */
            print_matrix(ddg_matrix, vectors_count, vectors_count);
        else if (strcmp(goal, "norm") == 0) {
            /* If the goal isn't to calculate the diagonal degree matrix, then it's to calculate the normalized similarity matrix. */
            normalized_similarity_matrix = calculate_normalized_similarity_matrix(sym_matrix, ddg_matrix);
            if(normalized_similarity_matrix == NULL)
                return 1;
            print_matrix(normalized_similarity_matrix, vectors_count, vectors_count);
            free_memory_of_matrix(normalized_similarity_matrix, vectors_count);
        }
        free_memory_of_matrix(ddg_matrix, vectors_count);
    }
    free_memory_of_matrix(sym_matrix, vectors_count);
    return 0;
}


/**
 * This function reads the file specified by the user to determine the number of vectors
 * and their dimensions, and then performs the specified goal operation (similarity matrix,
 * diagonal degree matrix, or normalized similarity matrix).
 * @param argc The number of command-line arguments.
 * @param argv An array of strings representing the command-line arguments.
 *             argv[1] should be the goal operation ("sym", "ddg", or "norm"),
 *             and argv[2] should be the input file name containing vector data.
 * @return 0 if the program executes successfully, or 1 if an error occurs.
 */
int main(int argc, char **argv){
    char *file_name, *goal;
    if (argc != 3) {
        printf(ERR_MSG);
        return 1;
    }
    file_name = argv[2];
    goal = argv[1];
    if(extract_data_from_file(file_name) != 0) {
        printf(ERR_MSG);
        return 1;
    }
    if(perform_goal(goal) != 0) {
        printf(ERR_MSG);
        return 1;
    }
    free_memory_of_matrix(data_vectors, vectors_count);
    return 0;
}

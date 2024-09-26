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

double calculate_squared_euclidean_distance(double *first_vector, double *second_vector)
{
    // calculating the squared Euclidean according to the formula
    int i;
    double sum_of_points = 0.0;
    for (i = 0 ; i < vector_length; i++)
    {
        sum_of_points += pow(first_vector[i] - second_vector[i], 2);
    }
    return sum_of_points;
}

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


void copy_vector_by_cord(double *copy_from, double *copy_to)
{
    // copy the array
    int i;
    for(i = 0; i < vector_length; i++)
    {
        copy_to[i] = copy_from[i];
    }
}


void free_memory_of_matrix(double **matrix, int numb_of_rows)
{
    // numb_of_rows- the number of rows that the matrix has or the number of rows created before failing
    // the method delete the rows' arrays one by one starting from the last row (numb_of_rows), and then delete the
    // array of the rows
    int i;
    for(i = 0; i < numb_of_rows; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

void free_memory_of_vectors_array(int vectors_counted)
{
    int i;
    for(i = 0; i < vectors_counted; i++)
    {
        free(data_vectors[i]);
    }
    free(data_vectors);
}


double **calculate_inverse_square_root(double **diagonal_matrix) {
    // receives the diagonal matrix which we want to find it's inverse square root, and create a new matrix such that
    // for every entry in the diagonal the value is (1/ sqrt(diagonal matrix received, in the same entry)) and for
    // every other entry is 0.
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

double **multiply_matrices(double **first_matrix, double **second_matrix, int rows_count, int columns_count, int common_count) {
    // the method multiply matrices according to the formula. the left matrix is first_matrix and the right one is second_matrix
    // moreover, the size of left one is (rows_count X common_count) and the size of the right one is (common_count X columns_count)
    // the matrix that the method return is the product of the multiplication, and it's size is (rows_count X columns_count)
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

double **transpose_matrix(double **matrix, int rows_count, int columns_count) {
    // switches the row and column indices of the matrix we want to transpose and put it in a new matrix
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

double **matrices_subtraction(double** first_matrix, double** second_matrix, int rows_count, int columns_count) {
    // subtract matrices by the formula: (A-B)(i,j) = a(i,j) - b(i,j)
    // and put it a new matrix, then return it
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

void print_matrix(double** matrix, int rows_count, int columns_count)
{
    // print the given matrix line by line
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
    // create array of rows
    sym_matrix = (double**) calloc(vectors_count, sizeof(double*));
    if (sym_matrix == NULL) {
        return NULL;
    }
    for (i = 0; i < vectors_count; i++){
        // create an array of row number i
        sym_matrix[i] = (double *)calloc( vectors_count, sizeof(double));
        if (sym_matrix[i] == NULL)
        {
            free_memory_of_matrix(sym_matrix, i);
            return NULL;
        }
        for (j = 0; j < vectors_count; j++){
            // calculate and insert to the entry the value of similarity matrix in the entry, the calculation is
            // according to the formula, using the "calculate_squared_euclidean_distance" method above
            if(i == j){
                sym_matrix[i][j] = 0.0;
            }
            else{
                sym_matrix[i][j] = exp(-0.5 * calculate_squared_euclidean_distance(data_vectors[i], data_vectors[j]));
            }
        }
    }
    return sym_matrix;
}

double **calculate_diagonal_degree_matrix(double **similarity_matrix)
{
    int i,j;
    double i_row_sum = 0.0;
    double **diagonal_degree_matrix;
    diagonal_degree_matrix = (double**) calloc(vectors_count, sizeof(double*));
    if (diagonal_degree_matrix == NULL) {
        return NULL; }
    for(i = 0; i < vectors_count; i++){
        diagonal_degree_matrix[i] = (double *)calloc(vectors_count, sizeof(double));
        i_row_sum = 0.0;
        if (diagonal_degree_matrix[i] == NULL)
        {
            free_memory_of_matrix(diagonal_degree_matrix, i);
            return NULL;
        }
        for(j = 0; j < vectors_count; j++) {
            i_row_sum += similarity_matrix[i][j];
        }
        diagonal_degree_matrix[i][i] = i_row_sum;
    }
    return diagonal_degree_matrix;
}

double **calculate_normalized_similarity_matrix(double** diagonal_degree_matrix , double** sym_matrix){
    double **inverse_square_root_matrix, **graph_Laplacian_left, **graph_Laplacian;
    // calculate the inverse square root matrix of the diagonal degree matrix that received
    inverse_square_root_matrix = calculate_inverse_square_root(diagonal_degree_matrix);
    if (inverse_square_root_matrix == NULL) {
        return NULL;
    }
    // multiply the inverse square root matrix of the diagonal degree matrix that received with the sym matrix received
    // such that the inverse square root matrix of the diagonal degree matrix is in the left and the sym matrix in the right
    graph_Laplacian_left = multiply_matrices(inverse_square_root_matrix, sym_matrix , vectors_count, vectors_count, vectors_count);
    if (graph_Laplacian_left == NULL) {
        free_memory_of_matrix(inverse_square_root_matrix, vectors_count);
        return NULL;
    }
    // multiply the product matrix calculated above with the inverse square root matrix of the diagonal degree matrix that received
    // such that the product matrix calculated above is in the left and the inverse square root matrix of the diagonal degree matrix in the right
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

double **calculate_H_Ht_H_matrix(double **decomposition_matrix_H ) {
    double **matrix_H_transposed, **matrix_H_H_t, **matrix_H_Ht_H;
    matrix_H_transposed = transpose_matrix(decomposition_matrix_H, K, vectors_count);
    if (matrix_H_transposed == NULL) {
        return NULL;
    }
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
            return NULL;
        }
        matrix_H_Ht_H = calculate_H_Ht_H_matrix(decomposition_matrix_H);
        if (matrix_H_Ht_H == NULL) {
            free_memory_of_matrix(decomposition_matrix_H, vectors_count);
            free_memory_of_matrix(matrix_WH, vectors_count);
        }
        updated_H = calculate_updated_H(decomposition_matrix_H, matrix_WH, matrix_H_Ht_H);
        if(updated_H == NULL) {
            free_memory_of_matrix(decomposition_matrix_H, vectors_count);
            return NULL;
        }
        subtraction_matrix = matrices_subtraction(updated_H, decomposition_matrix_H , vectors_count , K);
        if (subtraction_matrix == NULL) {
            free_memory_of_matrix(decomposition_matrix_H, vectors_count);
            free_memory_of_matrix(updated_H, vectors_count);
            return NULL;
        }
        frobenius_norm_value = calculate_frobenius_norm_squared(subtraction_matrix, vectors_count, K);
        free_memory_of_matrix(subtraction_matrix, vectors_count);
        free_memory_of_matrix(decomposition_matrix_H, vectors_count);
        decomposition_matrix_H = updated_H;
        if (frobenius_norm_value < EPSILON)
            break;
        iter_count++;
    }
    return decomposition_matrix_H;
}


int get_matrix_dimensions(FILE* vectors_file) {
    char line[MAX_LINE_LENGTH];
    vectors_count = 0;
    vector_length = 0;
    if (fgets(line, MAX_LINE_LENGTH, vectors_file)) {
        // read the first vector in order to find how many rows there are in each vector and does that by counting
        // these rows one by one
        char *token = strtok(line, ",");
        while (token != NULL) {
        vector_length++;
        token = strtok(NULL, ",");
        }
    }
    vectors_count++;
    while(fgets(line, MAX_LINE_LENGTH, vectors_file)) {
        // counting how many vectors there are in the given input starting from the second one (the first one counted above)
        vectors_count++;
    }

    // return to the beginning of the file in order to read the vectors afterward
    rewind(vectors_file);
    return 0;
}

int read_vectors(FILE* vectors_file) {
    int j, i;
    char line[MAX_LINE_LENGTH];
    char *token;
    i = 0;
    // create the vectors array
    data_vectors = (double**) calloc(vectors_count, sizeof(double*));
    if(data_vectors == NULL) {
        return 1;
    }
    while(fgets(line, MAX_LINE_LENGTH, vectors_file)) {
        // create the array of vector i
        data_vectors[i] = (double*) calloc(vector_length, sizeof(double));
        if(data_vectors[i] == NULL) {
            free_memory_of_matrix(data_vectors, i);
            return 1;
        }
        // read the values of the cells of vector i using the "," delimiter
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

// receive the input file name, find out how many vectors there are and how many rows each vector has, than read the file
// and insert the vectors to an array
int extract_data_from_file(FILE* file_name) {
    FILE* vectors_file;
    vectors_file = fopen(file_name, "r");
    if(vectors_file == NULL) {
        return 1;
    }
    get_matrix_dimensions(vectors_file);
    if(read_vectors(vectors_file) != 0) {
        return 1;
    }
    fclose(vectors_file);
    return 0;
}

// receive the goal given and perform it using the methods above
int perform_goal(char* goal) {
    double **sym_matrix, **ddg_matrix, **normalized_similarity_matrix;
    // calculate the similarity matrix of the vectors array matrix, the method will use it for any goal
    sym_matrix = calculate_similarity_matrix();
    if(sym_matrix == NULL) {
        return 1;
    }
    if(strcmp(goal, "sym") == 0) {
        // if the goal is to calculate the similarity matrix: print it and finish
        print_matrix(sym_matrix, vectors_count, vectors_count);
    }
    else {
        // if the goal isn't to calculate the similarity matrix, the method calculates the diagonal degree matrix
        // of  similarity matrix calculated above, the method will now use it for any goal
        ddg_matrix = calculate_diagonal_degree_matrix(sym_matrix);
        if(ddg_matrix == NULL) {
            free_memory_of_matrix(sym_matrix, vectors_count);
            return 1;
        }
        if(strcmp(goal, "ddg") == 0) {
            // if the goal is to calculate the diagonal degree matrix: print it and finish
            print_matrix(ddg_matrix, vectors_count, vectors_count);
        }
        else if (strcmp(goal, "norm") == 0) {
            // if the goal isn't to calculate the diagonal degree matrix, then it's to calculate the normalized similarity matrix
            normalized_similarity_matrix = calculate_normalized_similarity_matrix(sym_matrix, ddg_matrix);
            if(normalized_similarity_matrix == NULL) {
                return 1;
            }
            print_matrix(normalized_similarity_matrix, vectors_count, vectors_count);
            free_memory_of_matrix(normalized_similarity_matrix, vectors_count);
        }
        free_memory_of_matrix(ddg_matrix, vectors_count);
    }
    free_memory_of_matrix(sym_matrix, vectors_count);
    return 0;
}

// find out how many vectors there are and how many rows each vector has
// insert the vectors given in the input to an array using the methods "read_vectors"
//
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

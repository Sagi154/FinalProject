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


int vectors_count, vector_length, K;
double **data_vectors;
struct cluster *clusters = NULL;

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
    int i;
    for(i = 0; i < vector_length; i++)
    {
        copy_to[i] = copy_from[i];
    }
}

void free_vector_cords(struct cord *head_cord, int cords_counted)
{
    int i;
    struct cord *curr_cord;
    for ( i = 0; i < cords_counted; i++)
    {
        curr_cord = head_cord;
        head_cord = head_cord->next;
        free(curr_cord);
    }
}

void free_memory_of_lists(struct vector *head_vec, int vectors_counted)
{
    int i;
    struct vector *curr_vec;
    for ( i = 0; i < vectors_counted; i++)
    {
        curr_vec = head_vec;
        head_vec = head_vec->next;
        free_vector_cords(curr_vec->cords, vector_length);
        free(curr_vec);
    }
    if (head_vec != NULL)
    {
        free(head_vec);
    }
}

void free_memory_of_matrix(double **matrix, int numb_of_rows)
{
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

int create_vector_arr(struct vector *head_vec)
{
    int i, j;
    struct vector *tmp_vec;
    struct cord *tmp_cord;
    double *vector_i;
    tmp_vec = head_vec;
    data_vectors = (double**)calloc(vectors_count, sizeof(double*));
    if (data_vectors == NULL)
    {
        printf(ERR_MSG);
        return 1;
    }
    for (i = 0; i < vectors_count; i++)
    {
        tmp_cord = tmp_vec->cords;
        vector_i = (double *)calloc(vector_length, sizeof(double));
        if (vector_i == NULL)
        {
            printf(ERR_MSG);
            free_memory_of_vectors_array(i);
            return 1;
        }
        for (j = 0; j < vector_length; j++)
        {
            vector_i[j] = tmp_cord->value;
            tmp_cord = tmp_cord->next;
        }
        tmp_vec= tmp_vec->next;
        data_vectors[i] = vector_i;
    }
    return 0;
}

int is_double_integer(double value)
{
    return value == (int)value;
}

double **calculate_inverse_square_root(double **diagonal_matrix) {
    int i;
    double **inverse_square_root_matrix = (double **) calloc(vectors_count, sizeof(double*));
    if (inverse_square_root_matrix == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    for (i = 0; i < vectors_count; i++) {
        inverse_square_root_matrix[i] = (double *)calloc(vectors_count, sizeof(double));
        if (inverse_square_root_matrix[i] == NULL)
        {
            free_memory_of_matrix(inverse_square_root_matrix, i);
            printf(ERR_MSG);
            return NULL;
        }
        inverse_square_root_matrix[i][i] = 1/ sqrt(diagonal_matrix[i][i]);
    }
    return inverse_square_root_matrix;
}

double **multiply_matrices(double **first_matrix, double **second_matrix, int rows_count, int columns_count, int common_count) {
    int i,j,k;
    double **product_matrix = (double **) calloc(rows_count, sizeof(double *));
    if (product_matrix == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    for (i = 0; i < rows_count; i++) {
        product_matrix[i] = (double *)calloc(columns_count, sizeof(double));
        if (product_matrix[i] == NULL)
        {
            printf(ERR_MSG);
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
    int i,j;
    double **transpose_matrix = (double **) calloc(rows_count, sizeof(double *));
    if (transpose_matrix == NULL)
    {
        printf(ERR_MSG);
        return NULL;
    }
    for (i = 0; i < rows_count; i++) {
        transpose_matrix[i] = (double *)calloc(columns_count, sizeof(double));
        if (transpose_matrix[i] == NULL)
        {
            free_memory_of_matrix(transpose_matrix, i);
            printf(ERR_MSG);
            return NULL;
        }
        for (j = 0; j < columns_count; j++) {
            transpose_matrix[i][j] = matrix[j][i];
        }
    }
    return transpose_matrix;
}

double **matrices_subtraction(double** first_matrix, double** second_matrix, int rows_count, int columns_count) {
    int i,j;
    double **subtraction_matrix = (double **) calloc(rows_count, sizeof(double *));
    if (subtraction_matrix == NULL)
    {
        printf(ERR_MSG);
        return NULL;
    }
    for (i = 0; i < rows_count; i++) {
        subtraction_matrix[i] = (double *)calloc(columns_count, sizeof(double));
        if (subtraction_matrix[i] == NULL)
        {
            free_memory_of_matrix(subtraction_matrix, i);
            printf(ERR_MSG);
            return NULL;
        }
        for (j = 0; j < columns_count; j++) {
            subtraction_matrix[i][j] = first_matrix[i][j] - second_matrix[i][j];
        }
    }
    return subtraction_matrix;
}

double avg_W_entries(double ** normalized_similarity_matrix) {
    int i,j;
    double sum = 0;
    for (i = 0; i < vectors_count; i++) {
        for (j = 0; j < vectors_count; j++) {
            sum += normalized_similarity_matrix[i][j];
        }
    }
    return sum / (vectors_count * vectors_count);
}

double **initialize_H(double **normalized_similarity_matrix) {
    int i,j;
    double** decomposition_matrix;
    double m, bound;
    m = avg_W_entries(normalized_similarity_matrix);
    decomposition_matrix = (double**)calloc(vectors_count, sizeof(double*));
    if (decomposition_matrix == NULL){
        printf(ERR_MSG);
        return NULL; }
    bound = 2 * sqrt(m/K);
    for (i = 0; i < vectors_count; i++) {
        decomposition_matrix[i] = (double *)calloc( K, sizeof(double));
        if (decomposition_matrix[i] == NULL){
            printf(ERR_MSG);
            free_memory_of_matrix(decomposition_matrix, i);
            return NULL; }
        for (j = 0; j < K; j++) {
            decomposition_matrix[i][j] = ((double)rand() / RAND_MAX) * bound;
        }
    }
    return decomposition_matrix;
}

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
    if (sym_matrix == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    for (i = 0; i < vectors_count; i++){
        sym_matrix[i] = (double *)calloc( vectors_count, sizeof(double));
        if (sym_matrix[i] == NULL)
        {
            free_memory_of_matrix(sym_matrix, i);
            printf(ERR_MSG);
            return NULL;
        }
        for (j = 0; j < vectors_count; j++){
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
        printf(ERR_MSG);
        return NULL; }
    for(i = 0; i < vectors_count; i++){
        diagonal_degree_matrix[i] = (double *)calloc(vectors_count, sizeof(double));
        i_row_sum = 0.0;
        if (diagonal_degree_matrix[i] == NULL)
        {
            free_memory_of_matrix(diagonal_degree_matrix, i);
            printf(ERR_MSG);
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
    inverse_square_root_matrix = calculate_inverse_square_root(diagonal_degree_matrix);
    if (inverse_square_root_matrix == NULL) {
        return NULL;
    }

    graph_Laplacian_left = multiply_matrices(inverse_square_root_matrix, sym_matrix , vectors_count, vectors_count, vectors_count);
    if (graph_Laplacian_left == NULL) {
        free_memory_of_matrix(inverse_square_root_matrix, vectors_count);
        return NULL;
    }

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

double **calculate_final_decomposition_matrix_symnmf(double **decomposition_matrix_H, double **normalized_similarity_matrix) {
    int iter_count;
    double frobenius_norm_value;
    double **matrix_WH, **matrix_H_transposed, **matrix_H_H_t, **matrix_H_Ht_H, **subtraction_matrix;
    double **updated_H;
    iter_count = 1;
    while (iter_count <= MAX_ITER) {
        int i,j;
        matrix_WH = multiply_matrices(normalized_similarity_matrix, decomposition_matrix_H, vectors_count, K, vectors_count);
        if (matrix_WH == NULL) {
            free_memory_of_matrix(decomposition_matrix_H, vectors_count);
            return NULL;
        }
        matrix_H_transposed = transpose_matrix(decomposition_matrix_H, K, vectors_count);
        if (matrix_H_transposed == NULL) {
            free_memory_of_matrix(decomposition_matrix_H, vectors_count);
            free_memory_of_matrix(matrix_WH, vectors_count);
            return NULL;
        }
        matrix_H_H_t = multiply_matrices(decomposition_matrix_H, matrix_H_transposed, vectors_count, vectors_count, K);
        if (matrix_H_H_t == NULL) {
            free_memory_of_matrix(decomposition_matrix_H, vectors_count);
            free_memory_of_matrix(matrix_WH, vectors_count);
            free_memory_of_matrix(matrix_H_transposed, K);
            return NULL;
        }
        matrix_H_Ht_H = multiply_matrices(matrix_H_H_t, decomposition_matrix_H, vectors_count, K, vectors_count);
        if (matrix_H_Ht_H == NULL) {
            free_memory_of_matrix(decomposition_matrix_H, vectors_count);
            free_memory_of_matrix(matrix_WH, vectors_count);
            free_memory_of_matrix(matrix_H_transposed, K);
            free_memory_of_matrix(matrix_H_H_t, vectors_count);
            return NULL;
        }
        free_memory_of_matrix(matrix_H_transposed, K);
        free_memory_of_matrix(matrix_H_H_t, vectors_count);
        updated_H = (double**)calloc(vectors_count, sizeof(double*));
        if (updated_H == NULL) {
            printf(ERR_MSG);
            free_memory_of_matrix(decomposition_matrix_H, vectors_count);
            free_memory_of_matrix(matrix_WH, vectors_count);
            return NULL;
        }
        for ( i = 0; i < vectors_count; i++) {
            updated_H[i] = (double *)calloc(K, sizeof(double));
            if (updated_H[i] == NULL) {
                printf(ERR_MSG);
                free_memory_of_matrix(decomposition_matrix_H, vectors_count);
                free_memory_of_matrix(matrix_WH, vectors_count);
                free_memory_of_matrix(updated_H, i);
                return NULL;
            }
            for ( j = 0; j < K; j++) {
                updated_H[i][j] = decomposition_matrix_H[i][j] * (1 - BETA + BETA * (matrix_WH[i][j] / matrix_H_Ht_H[i][j]));
            }
        }
        free_memory_of_matrix(matrix_WH, vectors_count);
        free_memory_of_matrix(matrix_H_Ht_H, vectors_count);
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


int get_matrix_dimensions(char* file_name) {
    char c;
    int flag;
    FILE* vectors_file;
    vectors_count = 0;
    vector_length = 0;
    flag = 0;
    vectors_file = fopen(file_name, "r");
    if(vectors_file == NULL) {
        printf(ERR_MSG);
        return 1;
    }

    while((c = fgetc(vectors_file)) != EOF) {
        if(c == '\n') {
            if(flag == 0) {
                flag = 1;
            }
            vectors_count++;
        }
        if(flag == 0) {
            vector_length++;
        }
    }
    rewind(vectors_file);
    fclose(vectors_file);
    return 0;
}

int read_vectors(char* file_name) {
    int i,j;
    FILE* vectors_file;
    vectors_file = fopen(file_name, "r");
    if(vectors_file == NULL) {
        printf(ERR_MSG);
        return 1;
    }
    data_vectors = (double**) calloc(vectors_count, sizeof(double*));
    if(data_vectors == NULL) {
        printf(ERR_MSG);
        return 1;
    }
    for(i = 0; i < vectors_count; i++) {
        data_vectors[i] = (double*) calloc(vector_length, sizeof(double));
        if(data_vectors[i] == NULL) {
            printf(ERR_MSG);
            free_memory_of_matrix(data_vectors, i);
            return 1;
        }
        for(j = 0; j < vector_length; j++) {
            fscanf(vectors_file, "%lf", &data_vectors[i][j]);
        }
        fclose(vectors_file);
    }
    return 0;
}

int main(int argc, char **argv){
    char *file_name, *goal;
    double **sym_matrix, **ddg_matrix, **normalized_similarity_matrix;
    if (argc != 3) {
        printf(ERR_MSG);
        return 1;
    }
    file_name = argv[2];
    goal = argv[1];
    if(get_matrix_dimensions(file_name) == 1) {
        return 1;
    }
    if(read_vectors(file_name) == 1) {
        return 1;
    }
    printf("About to print data vectors: \n");
    printf("Vectors count: %d, Vector length: %d \n", vectors_count, vector_length);
    print_matrix(data_vectors, vectors_count, vector_length);
    sym_matrix = calculate_similarity_matrix();
    if(sym_matrix == NULL) {
        return 1;
    }
    if(strcmp(goal, "sym") == 0) {
        print_matrix(sym_matrix, vectors_count, vectors_count);
    }
    else {
        ddg_matrix = calculate_diagonal_degree_matrix(sym_matrix);
        if(ddg_matrix == NULL) {
            free_memory_of_matrix(sym_matrix, vectors_count);
            return 1;
        }
        if(strcmp(goal, "ddg") == 0) {
            print_matrix(ddg_matrix, vectors_count, vectors_count);
        }
        else if (strcmp(goal, "norm") == 0) {
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
    free_memory_of_matrix(data_vectors, vectors_count);
    return 0;
}

#include <symnmf.h>
#include <math.h>
#include <stdlib.h>

#define ERR_MSG "An Error Has Occurred"
#define DBL_MAX 1.7976931348623157e+308
#define BETA 0.5
#define EPSILON 0.0001
#define MAX_ITER 300


int vectors_count, vector_length, K, failure;
double **vector_array;
struct cluster *clusters = NULL;

double calculate_squared_euclidean_distance(double *first_vector, double *second_vector)
{
    double sum_of_points;
    int i;
    sum_of_points = 0.0;
    for (i = 0 ; i < vector_length; i++)
    {
        sum_of_points += pow(first_vector[i]- second_vector[i], 2);
    }
    return sum_of_points;
}

double calculate_frobenius_norm_squared(double **matrix, int rows_count, int columns_count)
{
    int i,j;
    double sum = 0.0;
    for (i = 0; i < rows_count; i++) {
        for (j = 0; j < columns_count; j++) {
            sum += matrix[i][j] * matrix[i][j];
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

void free_memory_of_matrix(double **matrix, int index)
{
    int i;
    for(i = 0; i < index; i++)
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
        free(vector_array[i]);
    }
    free(vector_array);
}

int create_vector_arr(struct vector *head_vec)
{
    int i, j;
    struct vector *tmp_vec;
    struct cord *tmp_cord;
    double *vector_i;
    tmp_vec = head_vec;
    vector_array = (double**)calloc(sizeof(double*), vectors_count);
    if (vector_array == NULL)
    {
        printf(ERR_MSG);
        return 1;
    }
    for (i = 0; i < vectors_count; i++)
    {
        tmp_cord = tmp_vec->cords;
        vector_i = (double *)calloc(sizeof(double), vector_length);
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
        vector_array[i] = vector_i;
    }
    return 0;
}

int is_double_integer(double value)
{
    return value == (int)value;
}

double **calculate_inverse_square_root(double **diagonal_matrix) {
    int i;
    double **inverse_square_root_matrix = (double **) calloc(sizeof(double *), vector_length);
    if (inverse_square_root_matrix == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    for (i = 0; i < vectors_count; i++) {
        inverse_square_root_matrix[i] = (double *)calloc(sizeof(double), vector_length);
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
    double **product_matrix = (double **) calloc(sizeof(double *), rows_count);
    if (product_matrix == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    for (i = 0; i < rows_count; i++) {
        product_matrix[i] = (double *)calloc(sizeof(double), columns_count);
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
    double **transpose_matrix = (double **) calloc(sizeof(double *), rows_count);
    if (transpose_matrix == NULL)
    {
        printf(ERR_MSG);
        return NULL;
    }
    for (i = 0; i < rows_count; i++) {
        transpose_matrix[i] = (double *)calloc(sizeof(double), columns_count);
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
    double **subtraction_matrix = (double **) calloc(sizeof(double *), rows_count);
    if (subtraction_matrix == NULL)
    {
        printf(ERR_MSG);
        return NULL;
    }
    for (i = 0; i < rows_count; i++) {
        subtraction_matrix[i] = (double *)calloc(sizeof(double), columns_count);
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
    for (int i = 0; i < vector_length; i++) {
        for (int j = 0; j < vector_length; j++) {
            sum += normalized_similarity_matrix[i][j];
        }
    }
    return sum / pow(vector_length, 2);
}

void initialize_H(double **decomposition_matrix, double m, int rows_count, int columns_count) {
    int i,j;
    double range = 2 * sqrt(m/K);
    for (i = 0; i < rows_count; i++) {
        for (j = 0; j < columns_count; j++) {
            decomposition_matrix[i][j] = ((double)rand() / RAND_MAX) * range;
        }
    }
}

double **sym(){
/*
Calculate and output the similarity matrix as described in 1.1.
 */
    int i,j;
    double **sym_matrix = (double**) calloc(sizeof(double*), vectors_count);
    if (sym_matrix == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    for (i = 0; i < vectors_count; i++){
        sym_matrix[i] = (double *)calloc(sizeof(double), vectors_count);
        if (sym_matrix[i] == NULL)
        {
            printf(ERR_MSG);
            return NULL;
        }
        for (j = 0; j < vectors_count; j++){
            if(i == j){
                sym_matrix[i][j] = 0;
            }
            else{
                sym_matrix[i][j] = exp(-0.5 * calculate_squared_euclidean_distance(vector_array[i], vector_array[j])) ;
            }
        }
    }
    return sym_matrix;
}

double **ddg(double ** similarity_matrix)
{
    int i,j;
    double i_row_sum = 0.0;
    double **diagonal_degree_matrix = (double**) calloc(vectors_count, sizeof(double*));
    for(i = 0; i < vectors_count; i++){
        diagonal_degree_matrix[i] = (double *)calloc(sizeof(double), vectors_count);
        for(j = 0; j < vectors_count; j++) {
            i_row_sum += similarity_matrix[i][j];
        }
        diagonal_degree_matrix[i][i] = i_row_sum;
    }
}

double **norm(double** diagonal_degree_matrix , double** sym_matrix){
    /*
    Calculate and output the normalized similarity matrix as described in 1.3
    */
    double** inverse_square_root_matrix = (double**) calloc(sizeof(double*), vectors_count);
    if (inverse_square_root_matrix == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    inverse_square_root_matrix = calculate_inverse_square_root(diagonal_degree_matrix);

    double** graph_Laplacian = (double**) calloc(sizeof(double*), vectors_count);
    if (graph_Laplacian == NULL) {
        printf(ERR_MSG);
        free_memory_of_matrix(inverse_square_root_matrix, vectors_count);
        return NULL;
    }

    graph_Laplacian = multiply_matrices(inverse_square_root_matrix, sym_matrix , vectors_count, vectors_count, vectors_count);
    if (graph_Laplacian == NULL) {
        free_memory_of_matrix(inverse_square_root_matrix, vectors_count);
        return NULL;
    }

    graph_Laplacian = multiply_matrices(graph_Laplacian, inverse_square_root_matrix, vectors_count, vectors_count, vectors_count);
    if (graph_Laplacian == NULL) {
        free_memory_of_matrix(inverse_square_root_matrix, vectors_count);
        free_memory_of_matrix(graph_Laplacian, vectors_count);
        return NULL;
    }

    free_memory_of_matrix(inverse_square_root_matrix, vectors_count);
    return graph_Laplacian;
}

double **symnmf(double **normalized_similarity_matrix) {
    /*
    Perform full the symNMF as described in 1 and output H.
    */
    int iter_count;
    double frobenius_norm_value, m;
    double **matrix_WH, **matrix_H_transposed, **matrix_H_H_t, **matrix_H_Ht_H, **subtraction_matrix;
    double **decomposition_matrix, **updated_H;
    iter_count = 0;
    m = avg_W_entries(normalized_similarity_matrix);
    decomposition_matrix = (double**)calloc(vectors_count, sizeof(double*));
    if (decomposition_matrix == NULL)
    {
        printf(ERR_MSG);
        return NULL;
    }
    initialize_H(decomposition_matrix, m, vectors_count, K);
    while (iter_count <= MAX_ITER) {
        int i,j;
        iter_count++;
        matrix_WH = multiply_matrices(normalized_similarity_matrix, decomposition_matrix, vectors_count, K, vectors_count);
        if (matrix_WH == NULL) {
            printf(ERR_MSG);
            free_memory_of_matrix(decomposition_matrix, vectors_count);
            return NULL;
        }
        matrix_H_transposed = transpose_matrix(decomposition_matrix, vectors_count, K);
        if (matrix_H_transposed == NULL) {
            printf(ERR_MSG);
            free_memory_of_matrix(decomposition_matrix, vectors_count);
            free_memory_of_matrix(matrix_WH, vectors_count);
            return NULL;
        }
        matrix_H_H_t = multiply_matrices(decomposition_matrix, matrix_H_transposed, vectors_count, vectors_count, K);
        if (matrix_H_H_t == NULL) {
            printf(ERR_MSG);
            free_memory_of_matrix(decomposition_matrix, vectors_count);
            free_memory_of_matrix(matrix_WH, vectors_count);
            free_memory_of_matrix(matrix_H_transposed, vectors_count);
            return NULL;
        }
        matrix_H_Ht_H = multiply_matrices(matrix_H_H_t, decomposition_matrix, vectors_count, K, vectors_count);
        if (matrix_H_Ht_H == NULL) {
            printf(ERR_MSG);
            free_memory_of_matrix(decomposition_matrix, vectors_count);
            free_memory_of_matrix(matrix_WH, vectors_count);
            free_memory_of_matrix(matrix_H_transposed, vectors_count);
            free_memory_of_matrix(matrix_H_H_t, vectors_count);
            return NULL;
        }
        free_memory_of_matrix(matrix_H_transposed, vectors_count);
        free_memory_of_matrix(matrix_H_H_t, vectors_count);
        updated_H = (double**)calloc(vectors_count, sizeof(double*));
        if (updated_H == NULL) {
            printf(ERR_MSG);
            free_memory_of_matrix(decomposition_matrix, vectors_count);
            free_memory_of_matrix(matrix_WH, vectors_count);
            return NULL;
        }
        for ( i = 0; i < vectors_count; i++) {
            updated_H[i] = (double *)calloc(sizeof(double), K);
            if (updated_H[i] == NULL) {
                printf(ERR_MSG);
                free_memory_of_matrix(decomposition_matrix, vectors_count);
                free_memory_of_matrix(matrix_WH, vectors_count);
                free_memory_of_matrix(updated_H, i);
                return NULL;
            }
            for ( j = 0; j < K; j++) {
                updated_H[i][j] = decomposition_matrix[i][j] * (1 - BETA + BETA * (matrix_WH[i][j] / matrix_H_Ht_H[i][j]));
            }
        }
        free_memory_of_matrix(matrix_WH, vectors_count);
        free_memory_of_matrix(matrix_H_Ht_H, vectors_count);
        subtraction_matrix = matrices_subtraction(updated_H, decomposition_matrix , vectors_count , K);
        if (subtraction_matrix == NULL) {
            printf(ERR_MSG);
            free_memory_of_matrix(decomposition_matrix, vectors_count);
            free_memory_of_matrix(updated_H, vectors_count);
            return NULL;
        }
        frobenius_norm_value = calculate_frobenius_norm_squared(subtraction_matrix, vectors_count, K);
        free_memory_of_matrix(subtraction_matrix, vectors_count);
        free_memory_of_matrix(decomposition_matrix, vectors_count);
        decomposition_matrix = updated_H;
        if (frobenius_norm_value < EPSILON)
            break;
        }
    return decomposition_matrix;
}

int main(int argc, char **argv)
{
    struct vector *head_vec, *curr_vec;
    struct cord *head_cord, *curr_cord;
    double n, k_arg, iter_limit_arg;
    char c, *endptr1, *endptr2;
    int flag, cord_count;
    vectors_count = 0;
    cord_count = 0;
    vector_length = 1;
    flag = 0;
    int iter_limit = 200;
    if (argc == 2)
    {
        k_arg = strtod(argv[1], &endptr1);
    }
    else if (argc == 3)
    {
        k_arg = strtod(argv[1], &endptr1);
        iter_limit_arg = strtod(argv[2], &endptr2);
    }

    if (*endptr1 != '\0') {
        printf("Invalid number of clusters!\n");
        return 1;
    }
    if (!is_double_integer(k_arg))
    {
        printf("Invalid number of clusters!\n");
        return 1;
    }
    K = (int)k_arg;
    if (!(1 < K ))
    {
        printf("Invalid number of clusters!\n");
        return 1;
    }
    if (argc == 3)
    {
        if (*endptr2 != '\0')
        {
            printf("Invalid maximum iteration!\n");
            return 1;
        }
        if (!is_double_integer(iter_limit_arg))
        {
            printf("Invalid maximum iteration!\n");
            return 1;
        }
        iter_limit = (int)iter_limit_arg;
        if ( !(1 < iter_limit && iter_limit < 1000))
        {
            printf("Invalid maximum iteration!\n");
            return 1;
        }
    }

    head_cord = calloc(sizeof(struct cord), 1);
    if (head_cord == NULL)
    {
        printf(ERR_MSG);
        return 1;
    }
    curr_cord = head_cord;
    curr_cord->next = NULL;
    cord_count++;
    head_vec = calloc(sizeof(struct vector), 1);
    if (head_vec == NULL)
    {
        printf(ERR_MSG);
        free(head_cord);
        return 1;
    }
    curr_vec = head_vec;
    curr_vec->next = NULL;
    while (scanf("%lf%c", &n, &c) == 2)
    {
        if (c == '\n')
        {
            cord_count = 0;
            flag = 1;
            vectors_count++;
            curr_cord->value = n;
            curr_vec->cords = head_cord;
            curr_vec->next = calloc(sizeof(struct vector), 1);
            if (curr_vec->next == NULL)
            {
                printf(ERR_MSG);
                free_memory_of_lists(head_vec, vectors_count);
                return 1;
            }
            curr_vec = curr_vec->next;
            curr_vec->next = NULL;
            head_cord = calloc(sizeof(struct cord), 1);
            if (head_cord == NULL)
            {
                printf(ERR_MSG);
                free_memory_of_lists(head_vec, vectors_count);
                return 1;
            }
            cord_count++;
            curr_cord = head_cord;
            curr_cord->next = NULL;
            continue;
        }
        if (!flag)
        {
        vector_length++;
        }
        curr_cord->value = n;
        curr_cord->next = calloc(sizeof(struct cord), 1);
        if (curr_cord->next == NULL)
        {
            printf(ERR_MSG);
            fprintf(stderr, ERR_MSG);
            free_vector_cords(head_cord, cord_count);
            free_memory_of_lists(head_vec, vectors_count);
            return 1;
        }
        cord_count++;
        curr_cord = curr_cord->next;
        curr_cord->next = NULL;
    }
    if (!(K < vectors_count))
    {
        printf("Invalid number of clusters!\n");
        free_memory_of_lists(head_vec, vectors_count);
        free(head_cord);
        return 1;
    }
    failure = create_vector_arr(head_vec);
    if (failure)
    {
        free_memory_of_lists(head_vec, vectors_count);
        free(head_cord);
        return 1;
    }
    free_memory_of_lists(head_vec, vectors_count);
    free(head_cord);
}

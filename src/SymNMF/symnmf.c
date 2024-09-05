#include <symnmf.h>
#include <math.h>
#include <stdlib.h>

#define ERR_MSG "An Error Has Occurred"
#define DBL_MAX 1.7976931348623157e+308

int vectors_count, vector_length, iter_limit, K, failure;
double** vector_array;
struct cluster* clusters = NULL;

double calculate_squared_euclidean_distance(double* first_vector, double* second_vector)
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

double calculate_frobenius_form(){
}

void copy_vector_by_cord(double* copy_from, double* copy_to)
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

void free_memory_of_vectors_array(int vectors_counted)
{
    int i;
    for(i = 0; i < vectors_counted; i++)
    {
        free(vector_array[i]);
    }
    free(vector_array);
}

void free_memory_of_matrix(double** matrix, int index)
{
    int i;
    for(i = 0; i < index; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

double** calculate_inverse_square_root(double** diagonal_matrix) {
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

double** multiply_matrices(double** first_matrix, double** second_matrix) {
    int i,j,k;
    double **product_matrix = (double **) calloc(sizeof(double *), vector_length);
    if (product_matrix == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    for (i = 0; i < vectors_count; i++) {
        product_matrix[i] = (double *)calloc(sizeof(double), vector_length);
        if (product_matrix[i] == NULL)
        {
            printf(ERR_MSG);
            free_memory_of_matrix(product_matrix, i);
            return NULL;
        }
        for (j = 0; j < vector_length; j++) {
            for (k = 0; k < vector_length; k++) {
                product_matrix[i][j] += first_matrix[i][k] * second_matrix[k][j];
            }
        }
    }
    return product_matrix;

}

int create_vector_arr(struct vector* head_vec)
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

double **symnmf();
/*
Perform full the symNMF as described in 1 and output H.
 */
double **sym(){
        int i,j;
        double** sym_matrix = (double**) calloc(sizeof(double*), vector_length);
        if (sym_matrix == NULL) {
            printf(ERR_MSG);
            return NULL;
        }
        for (i = 0; i < vectors_count; i++){
            sym_matrix[i] = (double *)calloc(sizeof(double), vector_length);
            if (sym_matrix[i] == NULL)
            {
                printf(ERR_MSG);
                return NULL;
            }
          for (j = 0; j < vector_length; j++){
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
/*
Calculate and output the similarity matrix as described in 1.1.
 */
double **ddg();
/*
Calculate and output the Diagonal Degree Matrix as described in 1.2
 */
double **norm(double** diagonal_degree_matrix , double** sym_matrix) {
    double** inverse_square_root_matrix = (double**) calloc(sizeof(double*), vector_length);
    if (inverse_square_root_matrix == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    inverse_square_root_matrix = calculate_inverse_square_root(diagonal_degree_matrix);

    double** graph_Laplacian = (double**) calloc(sizeof(double*), vector_length);
    if (graph_Laplacian == NULL) {
        printf(ERR_MSG);
        free_memory_of_matrix(inverse_square_root_matrix, vectors_count);
        return NULL;
    }

    graph_Laplacian = multiply_matrices(inverse_square_root_matrix, sym_matrix);
    if (graph_Laplacian == NULL) {
        free_memory_of_matrix(inverse_square_root_matrix, vectors_count);
        return NULL;
    }

    graph_Laplacian = multiply_matrices(graph_Laplacian, inverse_square_root_matrix);
    if (graph_Laplacian == NULL) {
        free_memory_of_matrix(inverse_square_root_matrix, vectors_count);
        free_memory_of_matrix(graph_Laplacian, vectors_count);
        return NULL;
    }

    free_memory_of_matrix(inverse_square_root_matrix, vectors_count);
    return graph_Laplacian;
}
/*
Calculate and output the normalized similarity matrix as described in 1.3
 */



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
    iter_limit = 200;
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

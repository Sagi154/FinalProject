#define PY_SSIZE_T_CLEAN
#include <Python.h>

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "symnmf.h"

# define ERR_MSG "An Error Has Occurred\n"
# define DBL_MAX 1.7976931348623157e+308

double **transform_python_matrix_to_c_matrix(PyObject *lst_matrix, int rows_count, int columns_count)
{
    int i, j;
    double value;
    PyObject* item;
    PyObject* vect;
    double **c_matrix;
    c_matrix = (double**)calloc(rows_count, sizeof(double*));
    if (c_matrix == NULL)
    {
        printf(ERR_MSG);
        return NULL;
    }
    for (i = 0; i < vectors_count; i++)
    {
        c_matrix[i] = (double *)calloc(columns_count, sizeof(double));
        if (c_matrix[i] == NULL)
        {
            printf(ERR_MSG);
            free_memory_of_matrix(c_matrix, i);
            return NULL;
        }
        vect = PyList_GetItem(lst_matrix, i);
        for (j = 0; j < columns_count; j++)
        {
            item = PyList_GetItem(vect, j);
            value = PyFloat_AsDouble(item);
            c_matrix[i][j] = value;
        }
    }
    return c_matrix;
}

PyObject *transform_c_matrix_to_python_matrix(double **c_matrix, int rows_count, int columns_count) {
    int i, j;
    PyObject *python_matrix, *vector_i;
    PyObject *item_j;
    python_matrix = PyList_New(rows_count);
    if (python_matrix == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    for (i = 0; i < rows_count; i++)
    {
        vector_i =  PyList_New(columns_count);
        if (!PyList_Check(vector_i))
        {
            printf(ERR_MSG);
            return NULL;
        }
        for (j = 0; j < columns_count; j++)
        {
            item_j = Py_BuildValue("d", c_matrix[i][j]);
            if (!item_j)
            {
                printf(ERR_MSG);
                return NULL;
            }
            PyList_SetItem(vector_i, j, item_j);
        }
        PyList_SetItem(python_matrix, i, vector_i);
    }
    return python_matrix;
}



static PyObject* sym(PyObject *self, PyObject *args) {
    PyObject *lst_data_points;
    double **sym_matrix;
    PyObject* python_sym_matrix;
    if(!PyArg_ParseTuple(args, "Oi", &lst_data_points, &vectors_count))
    {
        printf(ERR_MSG);
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    data_vectors = transform_python_matrix_to_c_matrix(lst_data_points, vectors_count, vectors_count);
    if (data_vectors == NULL)
    {
        printf(ERR_MSG);
        return NULL;
    }
    sym_matrix = calculate_similarity_matrix();
    free_memory_of_matrix(data_vectors, vectors_count);
    if (sym_matrix == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    python_sym_matrix = transform_c_matrix_to_python_matrix(sym_matrix, vectors_count, vectors_count);
    free_memory_of_matrix(sym_matrix, vectors_count);
    return python_sym_matrix;
}

static PyObject* ddg(PyObject *self, PyObject *args) {
    PyObject *lst_data_points;
    double **sym_matrix, **ddg_matrix;
    PyObject* python_ddg_matrix;
    if(!PyArg_ParseTuple(args, "Oi", &lst_data_points, &vectors_count))
    {
        printf(ERR_MSG);
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    data_vectors = transform_python_matrix_to_c_matrix(lst_data_points, vectors_count, vectors_count);
    if (data_vectors == NULL)
    {
        printf(ERR_MSG);
        return NULL;
    }
    sym_matrix = calculate_similarity_matrix();
    free_memory_of_matrix(data_vectors, vectors_count);
    if (sym_matrix == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    ddg_matrix = calculate_diagonal_degree_matrix(sym_matrix);
    free_memory_of_matrix(sym_matrix, vectors_count);
    if (ddg_matrix == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    python_ddg_matrix = transform_c_matrix_to_python_matrix(ddg_matrix, vectors_count, vectors_count);
    free_memory_of_matrix(ddg_matrix, vectors_count);
    return python_ddg_matrix;
}

static PyObject* norm(PyObject *self, PyObject *args) {
    PyObject *lst_data_points;
    double **sym_matrix, **ddg_matrix, **norm_matrix;
    PyObject* python_norm_matrix;
    if(!PyArg_ParseTuple(args, "Oi", &lst_data_points, &vectors_count))
    {
        printf(ERR_MSG);
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    data_vectors = transform_python_matrix_to_c_matrix(lst_data_points, vectors_count, vectors_count);
    if (data_vectors == NULL)
    {
        printf(ERR_MSG);
        return NULL;
    }
    sym_matrix = calculate_similarity_matrix();
    free_memory_of_matrix(data_vectors, vectors_count);
    if (sym_matrix == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    ddg_matrix = calculate_diagonal_degree_matrix(sym_matrix);
    if (ddg_matrix == NULL) {
        printf(ERR_MSG);
        free_memory_of_matrix(sym_matrix, vectors_count);
        return NULL;
    }
    norm_matrix = calculate_normalized_similarity_matrix(ddg_matrix, sym_matrix);
    free_memory_of_matrix(sym_matrix, vectors_count);
    free_memory_of_matrix(ddg_matrix, vectors_count);
    if (norm_matrix == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    python_norm_matrix = transform_c_matrix_to_python_matrix(norm_matrix, vectors_count, vectors_count);
    free_memory_of_matrix(norm_matrix, vectors_count);
    return python_norm_matrix;
}

static PyObject* symnmf(PyObject *self, PyObject *args) {
    PyObject *lst_decomposition_matrix_H, *lst_norm_matrix;
    double **initial_decomposition_matrix_H, **norm_matrix;
    double **decomposition_matrix_H;
    PyObject* python_optimal_decomposition_matrix;
    if(!PyArg_ParseTuple(args, "iOOi", &K, &lst_decomposition_matrix_H, &lst_norm_matrix, &vectors_count))
    {
        printf(ERR_MSG);
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    initial_decomposition_matrix_H = transform_python_matrix_to_c_matrix(lst_decomposition_matrix_H, vectors_count, K);
    if (initial_decomposition_matrix_H == NULL)
    {
        printf(ERR_MSG);
        return NULL;
    }
    norm_matrix = transform_python_matrix_to_c_matrix(lst_norm_matrix, vectors_count, vectors_count);
    if (norm_matrix == NULL)
    {
        printf(ERR_MSG);
        free_memory_of_matrix(initial_decomposition_matrix_H, vectors_count);
        return NULL;
    }
    decomposition_matrix_H = calculate_final_decomposition_matrix_symnmf(initial_decomposition_matrix_H, norm_matrix);
    free_memory_of_matrix(initial_decomposition_matrix_H, vectors_count);
    free_memory_of_matrix(norm_matrix, vectors_count);
    if (decomposition_matrix_H == NULL) {
        printf(ERR_MSG);
        return NULL;
    }
    python_optimal_decomposition_matrix = transform_c_matrix_to_python_matrix(decomposition_matrix_H, vectors_count, K);
    free_memory_of_matrix(decomposition_matrix_H, vectors_count);
    return python_optimal_decomposition_matrix;
}

static PyMethodDef symnmfMethods[] = {
    {"sym", (PyCFunction) sym, METH_VARARGS, PyDoc_STR("Calculate and output the similarity matrix")},
    {"ddg", (PyCFunction) ddg, METH_VARARGS, PyDoc_STR("Calculate and output the Diagonal Degree Matrix")},
    {"norm", (PyCFunction) norm, METH_VARARGS, PyDoc_STR("Calculate and output the normalized similarity matrix")},
    {"symnmf", (PyCFunction) symnmf, METH_VARARGS, PyDoc_STR("Perform the full symNMF")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef symnmf_module = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule",
    NULL,
    -1,
    symnmfMethods
};

PyMODINIT_FUNC PyInit_symnmfmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmf_module);
    if (!m) {
        return NULL;
    }
    return m;
}

int main(int argc, char **argv)
{
    return 0;
}
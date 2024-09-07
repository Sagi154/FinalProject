#define PY_SSIZE_T_CLEAN
#include <Python.h>

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "symnmf.h"

# define ERR_MSG "An Error Has Occurred\n"
# define DBL_MAX 1.7976931348623157e+308

static PyObject* sym(PyObject *self, PyObject *args) {

}

static PyObject* ddg(PyObject *self, PyObject *args) {

}

static PyObject* norm(PyObject *self, PyObject *args) {

}

static PyObject* symnmf(PyObject *self, PyObject *args) {

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
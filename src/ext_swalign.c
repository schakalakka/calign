//
// Created by andreas on 1/23/16.
//


#include <Python.h>
#include "swalign.h"

static PyObject *swalign_fast_smith_waterman(PyObject *self, PyObject *args) {
    const char *a;
    const char *b;
    int *max_score = malloc(sizeof(int) * 3);
    seq_pair problem;

    if (!PyArg_ParseTuple(args, "ss", &a, &b)) {
        return NULL;
    }

    problem.a = a;
    problem.alen = strlen(problem.a);
    problem.b = b;
    problem.blen = strlen(problem.b);

    max_score = fast_smith_waterman(&problem);
    return Py_BuildValue("ssiii", problem.a, problem.b, max_score[0], max_score[1], max_score[2]);
}


static PyMethodDef pyswalign_methods[] = {
        //"PythonName"  C0function name,    argument presentation,  description
        {"fast_smith_waterman", swalign_fast_smith_waterman, METH_VARARGS, "fast smith waterman return score only"},
        {NULL, NULL,                                         0, NULL}  /*Sentinel*/
};

static struct PyModuleDef extswalign = {
        PyModuleDef_HEAD_INIT,
        "pyswalign", /* name of module */
        "",          /* module documentation, may be NULL */
        -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        pyswalign_methods
};

PyMODINIT_FUNC
PyInit_ext_swalign(void) {
    return PyModule_Create(&extswalign);
}
//
// Created by andreas on 1/23/16.
//


#include <Python.h>
#include "calign.h"

static PyObject *calign_align(PyObject *self, PyObject *args) {
    const char *a;
    const char *b;
    const char *align_type; //global, semiglobal=semi-global or local
    int *max_score = malloc(sizeof(int) * 3);
    seq_pair problem;

    if (!PyArg_ParseTuple(args, "sss", &a, &b, &align_type)) {
        return NULL;
    }

    problem.a = a;
    problem.alen = strlen(problem.a);
    problem.b = b;
    problem.blen = strlen(problem.b);

    max_score = alignment(&problem, align_type);
    return Py_BuildValue("ssiii", problem.a, problem.b, max_score[0], max_score[1], max_score[2]);
}


static PyMethodDef pycalign_methods[] = {
        //"PythonName"  C0function name,    argument presentation,  description
        {"align", calign_align, METH_VARARGS, "return score only"},
        {NULL, NULL,            0, NULL}  /*Sentinel*/
};

static struct PyModuleDef extcalign = {
        PyModuleDef_HEAD_INIT,
        "pycalign", /* name of module */
        "",          /* module documentation, may be NULL */
        -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        pycalign_methods
};

PyMODINIT_FUNC
PyInit_ext_calign(void) {
    return PyModule_Create(&extcalign);
}
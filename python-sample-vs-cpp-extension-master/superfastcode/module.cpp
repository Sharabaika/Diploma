#include <Python.h>
#include <Windows.h>


PyObject* test_fun_impl(PyObject*, PyObject* args) 
{
	PyObject* arg1;
	PyObject* arg2;

	PyArg_UnpackTuple(args, "ref", 1, 2, &arg1, &arg2);

	double n1 = PyFloat_AsDouble(arg1) + 1.0f;
	double n2 = PyFloat_AsDouble(arg2) + 1.0f;

	PyObject* r1 = PyFloat_FromDouble(n1);
	PyObject* r2 = PyFloat_FromDouble(n2);

	PyObject* tuple = PyTuple_New(2);
	PyTuple_SetItem(tuple, 0, r1);
	PyTuple_SetItem(tuple, 1, r2);

	return tuple;
}

static PyMethodDef superfastcode_methods[] = {
	{ "test_fun", (PyCFunction)test_fun_impl, METH_O, nullptr},

	// Terminate the array with an object containing nulls.
	{ nullptr, nullptr, 0, nullptr }
};

static PyModuleDef superfastcode_module = {
	PyModuleDef_HEAD_INIT,
	"superfastcode",                        // Module name to use with Python import statements
	"Provides some functions, but faster",  // Module description
	0,
	superfastcode_methods                   // Structure that defines the methods of the module
};

PyMODINIT_FUNC PyInit_superfastcode() {
	return PyModule_Create(&superfastcode_module);
}
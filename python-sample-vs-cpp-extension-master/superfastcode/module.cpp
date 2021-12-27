#include <Python.h>
#include <Windows.h>
#include <vector>
#include "Converters.h"

using std::vector;
 
PyObject* FluidSolver_Wrapper(PyObject* _, PyObject* args)
{
	PyObject* X;
	PyObject* Y;
	PyObject* Nodes;

	PyArg_UnpackTuple(args, "ref", 1, 3, &X, &Y, &Nodes);

	vector<double> X_arr = listTupleToVector_Double(X);
	vector<double> Y_arr = listTupleToVector_Double(Y);
	vector<vector<double>> Nodes_arr = listListToVector_Double(Nodes);

	for (size_t i = 0; i < X_arr.size(); i++)
	{
		X_arr[i] = X_arr[i] + 1.0f;
	}

	for (size_t i = 0; i < Nodes_arr.size(); i++)
	{
		Nodes_arr[i][0] = i;
	}

	X = vectorToList_Double(X_arr);
	Y = vectorToList_Double(Y_arr);


	Nodes = vectorVectorToTuple_Double(Nodes_arr);


	PyObject* res = PyTuple_New(3);
	PyTuple_SetItem(res, 0, X);
	PyTuple_SetItem(res, 1, Y);
	PyTuple_SetItem(res, 2, Nodes);

	return res;

}

PyObject* test_fun_impl(PyObject*, PyObject* args) 
{
	PyObject* arg1;
	PyObject* arg2;

	PyArg_UnpackTuple(args, "ref", 1, 2, &arg1, &arg2);

	double n1 = PyFloat_AsDouble(arg1) + 1.0f;
	double n2 = PyFloat_AsDouble(arg2) + 2.0f;

	PyObject* r1 = PyFloat_FromDouble(n1);
	PyObject* r2 = PyFloat_FromDouble(n2);

	PyObject* tuple = PyTuple_New(2);
	PyTuple_SetItem(tuple, 0, r1);
	PyTuple_SetItem(tuple, 1, r2);

	return tuple;
}

static PyMethodDef superfastcode_methods[] = {
	{ "test_fun", (PyCFunction)test_fun_impl, METH_O, nullptr},
	{ "SolveFluids", (PyCFunction)FluidSolver_Wrapper, METH_O, nullptr},

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
#include <Python.h>
#include <Windows.h>
#include <vector>
#include "Converters.h"
#include "Solver.h"

using std::vector;
 
PyObject* FluidSolver_Wrapper(PyObject* _, PyObject* args)
{
	PyObject *x, *y, *triangles, *segment_indecies, *trig_neighbours, *node_neighbours;
	PyObject *re, *vx;
	PyObject *qPsi, * qW, *max_error, *max_cycles;

	PyArg_UnpackTuple(args, "ref", 1, 12, &x, &y, &triangles, &segment_indecies, &trig_neighbours, &node_neighbours, &re, &vx, &qPsi, &qW, &max_error, &max_cycles);

	Arr X_arr = listTupleToVector_Double(x);
	Arr Y_arr = listTupleToVector_Double(y);
	JaggedArr Triangles_arr = listListToVector_Int(triangles);
	JaggedArr Segments = listListToVector_Int(segment_indecies);
	JaggedArr Trig_neighbours = listListToVector_Int(trig_neighbours);
	JaggedArr Node_neighbours = listListToVector_Int(node_neighbours);

	double Re = PyFloat_AsDouble(re);
	double Vx = PyFloat_AsDouble(vx);

	double QPsi = PyFloat_AsDouble(qPsi);
	double QW = PyFloat_AsDouble(qW);
	double Max_error = PyFloat_AsDouble(max_error);
	int Max_cycles = PyFloat_AsDouble(max_cycles);

	Arr PsiRes;
	Arr WRes;

	Arr Psi_Error;
	Arr W_Error;

	SolveFluid_Implementation(X_arr, Y_arr, Triangles_arr, Segments, Trig_neighbours, Node_neighbours, Re, Vx, QPsi, QW, Max_error, Max_cycles, PsiRes, WRes, Psi_Error, W_Error);

	PyObject* res = PyTuple_New(4);
	PyTuple_SetItem(res, 0, vectorToList_Double(PsiRes));
	PyTuple_SetItem(res, 1, vectorToList_Double(WRes));
	PyTuple_SetItem(res, 2, vectorToList_Double(Psi_Error));
	PyTuple_SetItem(res, 3, vectorToList_Double(W_Error));

	return res;
}

PyObject* SolveFast_Wrapper(PyObject* _, PyObject* args)
{
	PyObject* x, * y, * triangles, * segment_indices, * trig_neighbours, * node_neighbours, *triangle_indeces;
	PyObject* fluid_domain_nodes_indeces_array, * is_a_fluid_region_array, * is_a_wall_array;
	PyObject* pr, * ra, * ram, * chi0;

	PyObject* qPsi, * qW, * qT;
	PyObject* max_error, * max_cycles;

	PyObject* psi, * w, * t;

	PyObject* h_triangles;
	PyObject* dHdx, * dHdy;

	PyObject* normal_x, * normal_y;

	PyArg_UnpackTuple(args, "ref", 1, 27,
		&x, &y, &triangles, &segment_indices, &trig_neighbours, &node_neighbours, &triangle_indeces,
		&fluid_domain_nodes_indeces_array, &is_a_fluid_region_array, &is_a_wall_array,
		&pr, &ra, &ram, &chi0,
		&qPsi, &qW, &qT,
		&max_error, &max_cycles,
		&psi, &w, &t,
		&h_triangles,
		&dHdx, &dHdy,
		&normal_x, &normal_y);

	Arr X_arr = listTupleToVector_Double(x);
	Arr Y_arr = listTupleToVector_Double(y);
	JaggedArr Triangles_arr = listListToVector_Int(triangles);
	ArrInt Segments = listTupleToVector_Int(segment_indices);
	JaggedArr Trig_neighbours = listListToVector_Int(trig_neighbours);
	JaggedArr Node_neighbours = listListToVector_Int(node_neighbours);
	ArrInt Triangle_indeces = listTupleToVector_Int(triangle_indeces);

	ArrInt Fluid_domain_nodes_indeces_array = listTupleToVector_Int(fluid_domain_nodes_indeces_array);
	ArrInt Is_a_fluid_region_array = listTupleToVector_Int(is_a_fluid_region_array);
	ArrInt Is_a_wall_array = listTupleToVector_Int(is_a_wall_array);

	double Pr= PyFloat_AsDouble(pr);
	double Ra = PyFloat_AsDouble(ra);
	double Ram = PyFloat_AsDouble(ram);
	double Chi0 = PyFloat_AsDouble(chi0);

	double QPsi = PyFloat_AsDouble(qPsi);
	double QW = PyFloat_AsDouble(qW);
	double QT = PyFloat_AsDouble(qT);

	double Max_error = PyFloat_AsDouble(max_error);
	int Max_cycles = PyFloat_AsDouble(max_cycles);

	Arr Psi = listTupleToVector_Double(psi);
	Arr W = listTupleToVector_Double(w);
	Arr T = listTupleToVector_Double(t);

	Arr H_triangles = listTupleToVector_Double(h_triangles);
	Arr DHdx = listTupleToVector_Double(dHdx);
	Arr DHdy = listTupleToVector_Double(dHdy);

	Arr Normal_x = listTupleToVector_Double(normal_x);
	Arr Normal_y = listTupleToVector_Double(normal_y);

	SolveFast_Implementation(
		X_arr, Y_arr, Triangles_arr, Segments, Trig_neighbours, Node_neighbours, triangle_indeces,
		Re, Vx, QPsi, QW, Max_error, Max_cycles, PsiRes, WRes, Psi_Error, W_Error);

	PyObject* res = PyTuple_New(4);
	PyTuple_SetItem(res, 0, vectorToList_Double(PsiRes));
	PyTuple_SetItem(res, 1, vectorToList_Double(WRes));
	PyTuple_SetItem(res, 2, vectorToList_Double(Psi_Error));
	PyTuple_SetItem(res, 3, vectorToList_Double(W_Error));

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
	{ "SolveFast", (PyCFunction)SolveFast_Wrapper, METH_O, nullptr},
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
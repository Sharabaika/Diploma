#ifndef PYUTILS_H__
#define PYUTILS_H__

#include <Python.h>
#include <vector>
using namespace std;

PyObject* vectorToList_Double(const vector<double>& data);

PyObject* vectorToTuple_Double(const vector<float>& data);

PyObject* vectorVectorToTuple_Double(const vector< vector< double > >& data);

vector<double> listTupleToVector_Double(PyObject* incoming);

vector<int> listTupleToVector_Int(PyObject* incoming);

vector<vector<double>> listListToVector_Double(PyObject* in);

#endif
#include <pybind11/pybind11.h>

unsigned int fibinacci(const unsigned int n);

namespace py = pybind11;

PYBIND11_MODULE(pybind_11_example, mod) {
    mod.def("fibinacci_cpp", &fibinacci, "Recursive fibinacci algorithm.");
}

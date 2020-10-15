#pragma once

#include "graph_matching/matching_problem_input.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace LPMP {
    namespace py_helper {

    void add_assignments(LPMP::graph_matching_input& instance, const py::array_t<double> assignments);

    void add_quadratic_terms(LPMP::graph_matching_input& instance, const py::array_t<double> quadratic_terms, const py::array_t<int> left_edges, const py::array_t<int> right_edges);

    void construct_from_arrays(LPMP::graph_matching_input& instance, const py::array_t<double> assignments, const py::array_t<double> quadratic_terms, const py::array_t<int> left_edges, const py::array_t<int> right_edges);

    py::array_t<char> get_assignment_mask(LPMP::graph_matching_input::labeling labeling, const py::array_t<double> assignments);

    py::array_t<char> get_quadratic_terms_mask(LPMP::graph_matching_input::labeling labeling, const py::array_t<double> quadratic_terms, const py::array_t<int> left_edges, const py::array_t<int> right_edges);

    }
}

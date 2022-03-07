#include "multicut/multicut_instance.h"
#include "multicut/multicut_greedy_additive_edge_contraction.h"
#include "multicut/multicut_greedy_edge_fixation.h"
#include "multicut/multicut_balanced_edge_contraction.h"
#include "multicut/multicut_text_input.h"
#include <fstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

void construct_instance(LPMP::multicut_instance& instance, const Eigen::Array<size_t, Eigen::Dynamic, 2>& edge_indices, const Eigen::Array<double, Eigen::Dynamic, 1>& edge_costs)
{

    for(size_t e=0; e<edge_costs.rows(); ++e)
    {
        const size_t i = edge_indices(e, 0);
        const size_t j = edge_indices(e, 1);
        const double w = edge_costs(e);
        instance.add_edge(i, j, w);
    }
}

py::array_t<char> get_edge_mask(const LPMP::multicut_instance& instance, const LPMP::multicut_edge_labeling& edge_labels)
{
    char* edge_mask = new char[instance.no_edges()];
    assert(instance.no_edges() == edge_labels.size());
    for(size_t e=0; e<instance.no_edges(); ++e)
        edge_mask[e] = edge_labels[e];

    return py::array({pybind11::ssize_t(instance.no_edges())}, edge_mask); 
}

py::array_t<int> get_cc_ids_mask(const LPMP::multicut_instance& instance, const std::vector<int>& node_labeling)
{
    int* cc_mask = new int[instance.no_nodes()];
    assert(instance.no_nodes() == node_labeling.size());
    for(size_t i=0; i<instance.no_nodes(); ++i)
    {
        cc_mask[i] = node_labeling[i];
    }

    return py::array({pybind11::ssize_t(instance.no_nodes())}, cc_mask); 
} 

PYBIND11_MODULE(multicut_py, m) {
    m.doc() = "python binding for LPMP multicut primal algorithms";

    py::class_<LPMP::multicut_edge_labeling>(m, "multicut_edge_labeling")
        .def(py::init<>()); 

    py::class_<LPMP::multicut_instance>(m, "multicut_instance")
        .def(py::init<>())
        .def(py::init([](const Eigen::Array<size_t, Eigen::Dynamic, 2>& edge_indices, const Eigen::Array<double, Eigen::Dynamic, 1>& edge_costs) {
                    LPMP::multicut_instance instance;
                    construct_instance(instance, edge_indices, edge_costs);
                    return instance;
                    }))
        .def("result_mask", [](const LPMP::multicut_instance& instance, const LPMP::multicut_edge_labeling& edge_labeling, const std::vector<int>& node_labeling) {
                return std::make_pair(
                        get_edge_mask(instance, edge_labeling),
                        get_cc_ids_mask(instance, node_labeling)
                        ); 
                });
        m.def("read", [](LPMP::multicut_instance& instance, const std::string& filename)
                {
                instance = LPMP::multicut_text_input::parse_file(filename);
                });

        m.def("multicut_gaec", [](const LPMP::multicut_instance& instance) {
                return LPMP::greedy_additive_edge_contraction_impl(instance);
                });
        m.def("multicut_gef", [](const LPMP::multicut_instance& instance) {
                return LPMP::multicut_greedy_edge_fixation_impl(instance);
                });
        m.def("multicut_bec", [](const LPMP::multicut_instance& instance) {
                return LPMP::multicut_balanced_edge_contraction_impl(instance);
                });
}

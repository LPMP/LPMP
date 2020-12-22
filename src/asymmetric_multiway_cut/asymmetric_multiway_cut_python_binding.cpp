#include "asymmetric_multiway_cut/asymmetric_multiway_cut_instance.h"
#include "asymmetric_multiway_cut/asymmetric_multiway_cut_gaec.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

void construct_instance(LPMP::asymmetric_multiway_cut_instance& instance, const std::vector<std::tuple<size_t,size_t,double>>& edge_costs, const py::array_t<double> node_costs)
{

    for(size_t e=0; e<edge_costs.size(); ++e)
    {
        const auto ec =  edge_costs[e];
        const size_t i = std::get<0>(ec);
        const size_t j = std::get<1>(ec);
        const double w = std::get<2>(ec);
        instance.edge_costs.add_edge(i, j, w);
    }

    const auto _node_costs = node_costs.unchecked<2>();
    const size_t nr_nodes = _node_costs.shape(0);
    const size_t nr_labels = _node_costs.shape(1);
    std::vector<double> costs;
    for(size_t i=0; i<nr_nodes; ++i)
    {
        costs.clear();
        for(size_t l=0; l<nr_labels; ++l)
        {
            costs.push_back(_node_costs(i,l)); 
        } 
        instance.node_costs.push_back(costs.begin(), costs.end());
    } 
}

py::array_t<char> get_edge_mask(const LPMP::asymmetric_multiway_cut_instance& instance, const LPMP::asymmetric_multiway_cut_labeling& labeling)
{
    char* edge_mask = new char[instance.nr_edges()];
    assert(instance.edge_costs.no_edges() == labeling.edge_labels.size());
    for(size_t e=0; e<instance.nr_edges(); ++e)
        edge_mask[e] = labeling.edge_labels[e];

    return py::array({instance.nr_edges()}, edge_mask); 
} 

py::array_t<char> get_label_mask(const LPMP::asymmetric_multiway_cut_instance& instance, const LPMP::asymmetric_multiway_cut_labeling& labeling)
{
    char* label_mask = new char[instance.nr_nodes()*instance.nr_labels()];
    assert(instance.nr_nodes() == labeling.node_labels.size());
    for(size_t i=0; i<instance.nr_nodes(); ++i)
    {
        for(size_t l=0; l<instance.nr_labels(); ++l)
        {
            assert(labeling.node_labels[i] < instance.nr_labels());
            label_mask[i*instance.nr_labels() + l] = labeling.node_labels[i] == l;
        }
    }

    return py::array({instance.nr_nodes(), instance.nr_labels()}, label_mask); 
} 

py::array_t<int> get_cc_ids_mask(const LPMP::asymmetric_multiway_cut_instance& instance, const LPMP::asymmetric_multiway_cut_labeling& labeling)
{
    int* cc_mask = new int[instance.nr_nodes()];
    assert(instance.nr_nodes() == labeling.node_connected_components_ids.size());
    for(size_t i=0; i<instance.nr_nodes(); ++i)
    {
        cc_mask[i] = labeling.node_connected_components_ids[i];
    }

    return py::array({instance.nr_nodes()}, cc_mask); 
} 


PYBIND11_MODULE(asymmetric_multiway_cut_py, m) {
    m.doc() = "python binding for LPMP asymmetric multiway cut";

    py::class_<LPMP::asymmetric_multiway_cut_labeling>(m, "asymmetric_multiway_cut_labeling")
        .def(py::init<>()); 

    py::class_<LPMP::asymmetric_multiway_cut_instance>(m, "asymmetric_multiway_cut_instance")
        .def(py::init<>())
        .def(py::init([](const std::vector<std::tuple<size_t,size_t,double>>& edge_costs, const py::array_t<double> node_costs) {
                    LPMP::asymmetric_multiway_cut_instance instance;
                    construct_instance(instance, edge_costs, node_costs);
                    return instance;
                    }))
        .def("evaluate", &LPMP::asymmetric_multiway_cut_instance::evaluate)
        .def("result_mask", [](const LPMP::asymmetric_multiway_cut_instance& instance, const LPMP::asymmetric_multiway_cut_labeling& labeling) {
                return std::make_tuple(
                        get_edge_mask(instance, labeling),
                        get_label_mask(instance, labeling),
                        get_cc_ids_mask(instance, labeling)
                        ); 
                });

        m.def("asymmetric_multiway_cut_gaec", [](const LPMP::asymmetric_multiway_cut_instance& instance) {
                return LPMP::asymmetric_multiway_cut_gaec(instance);
                });
}

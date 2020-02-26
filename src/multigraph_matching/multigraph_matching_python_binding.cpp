#include "graph_matching/matching_problem_input.h"
#include "multigraph_matching/multigraph_matching.hxx"
#include "graph_matching/graph_matching_python_binding_helper.h"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

// there are two versions of the python binding: One with varying edges and with constant edge set per graph

// varying edges
void construct_from_arrays(LPMP::multigraph_matching_input& instance, const std::vector<py::array_t<double>> assignments, const std::vector<py::array_t<double>> quadratic_terms, const std::vector<py::array_t<int>> left_edges, const std::vector<py::array_t<int>> right_edges)
{
    const std::size_t nr_problems = assignments.size();

    if(nr_problems != quadratic_terms.size())
        throw std::runtime_error("#quadratic terms != #assignments");
    if(nr_problems != left_edges.size())
        throw std::runtime_error("#left edges != #assignments");
    if(nr_problems != right_edges.size())
        throw std::runtime_error("#right edges != #assignments");

    const std::size_t nr_graphs = 0.5 + std::sqrt(0.25+2*nr_problems);//(+1 + std::sqrt(-1+4*2*nr_problems) )/2;
    if((nr_graphs*(nr_graphs-1))/2 != nr_problems)
        throw std::runtime_error("there must be k*(k-1)/2 pairwise graph matching problems");

    std::size_t c=0;
    for(std::size_t i=0; i<nr_graphs; ++i) {
        for(std::size_t j=i+1; j<nr_graphs; ++j) {
            LPMP::graph_matching_input gm;
            LPMP::py_helper::add_assignments(gm, assignments[c]);
            instance.push_back({i, j, gm}); 
            ++c;
        }
    }
}

std::vector<py::array_t<char>> get_assignment_masks(LPMP::multigraph_matching_input::labeling labeling, const std::vector<py::array_t<double>> assignments)
{
    const std::size_t nr_problems = assignments.size();
    if(labeling.size() != assignments.size())
        throw std::runtime_error("#assignments != #labelings");

    std::vector<py::array_t<char>> masks;
    for(std::size_t c=0; c<nr_problems; ++c)
        masks.push_back(LPMP::py_helper::get_assignment_mask(labeling[c].labeling, assignments[c]));

    return masks; 
}

std::vector<py::array_t<char>> get_quadratic_terms_masks(LPMP::multigraph_matching_input::labeling labeling, const std::vector<py::array_t<double>> quadratic_terms, const std::vector<py::array_t<int>> left_edges, const std::vector<py::array_t<int>> right_edges)
{
    const std::size_t nr_problems = quadratic_terms.size();
    if(left_edges.size() != nr_problems)
        throw std::runtime_error("#left edges != #assignments");
    if(nr_problems != right_edges.size())
        throw std::runtime_error("#right edges != #assignments");

    std::vector<py::array_t<char>> masks;
    for(std::size_t c=0; c<nr_problems; ++c)
        masks.push_back(LPMP::py_helper::get_quadratic_terms_mask(labeling[c].labeling, quadratic_terms[c], left_edges[c], right_edges[c]));

    return masks; 
}

// constant edge set
void construct_from_arrays(LPMP::multigraph_matching_input& instance, const std::vector<py::array_t<double>> assignments, const std::vector<py::array_t<double>> quadratic_terms, const std::vector<py::array_t<int>> edges)
{
    const std::size_t nr_problems = assignments.size();

    if(nr_problems != quadratic_terms.size())
        throw std::runtime_error("#quadratic terms != #assignments");

    const std::size_t nr_graphs = 0.5 + std::sqrt(0.25+2*nr_problems);//(+1 + std::sqrt(-1+4*2*nr_problems) )/2;
    if((nr_graphs*(nr_graphs-1))/2 != nr_problems)
        throw std::runtime_error("there must be k*(k-1)/2 pairwise graph matching problems");
    if(nr_graphs != edges.size())
        throw std::runtime_error("#left edges != #graphs");

    std::size_t c=0;
    for(std::size_t i=0; i<nr_graphs; ++i) {
        for(std::size_t j=i+1; j<nr_graphs; ++j) {
            LPMP::graph_matching_input gm;
            LPMP::py_helper::add_assignments(gm, assignments[c]);
            LPMP::py_helper::add_quadratic_terms(gm, quadratic_terms[c], edges[i], edges[j]);
            instance.push_back({i, j, gm}); 
            ++c;
        }
    }
}

std::vector<py::array_t<char>> get_quadratic_terms_mask(LPMP::multigraph_matching_input::labeling labeling, const std::vector<py::array_t<double>> quadratic_terms, const std::vector<py::array_t<int>> edges)
{
    const std::size_t nr_problems = quadratic_terms.size();
    if(labeling.size() != quadratic_terms.size())
        throw std::runtime_error("#labelings != #quadratic_terms");

    const std::size_t nr_graphs = 0.5 + std::sqrt(0.25+2*nr_problems);//(+1 + std::sqrt(-1+4*2*nr_problems) )/2;
    if((nr_graphs*(nr_graphs-1))/2 != nr_problems)
        throw std::runtime_error("there must be k*(k-1)/2 pairwise graph matching problems");
    if(nr_graphs != edges.size())
        throw std::runtime_error("#left edges != #graphs");

    std::vector<py::array_t<char>> masks;
    std::size_t c=0;
    for(std::size_t i=0; i<nr_graphs; ++i) {
        for(std::size_t j=i+1; j<nr_graphs; ++j) {
            if(labeling[c].left_graph_no != i || labeling[c].right_graph_no != j) 
                throw std::runtime_error("labelings not in lexicographical order");
            masks.push_back(LPMP::py_helper::get_quadratic_terms_mask(labeling[c].labeling, quadratic_terms[c], edges[i], edges[j]));
            ++c;
        }
    }

    return masks; 
}


PYBIND11_MODULE(multigraph_matching_py, m) {

    m.doc() = "python binding for LPMP multigraph matching";

    py::class_<LPMP::multigraph_matching_input::labeling>(m, "multigraph_matching_labeling")
        .def(py::init<>())
        .def("result_masks", [](const LPMP::multigraph_matching_input::labeling& l, const std::vector<py::array_t<double>> assignments, const std::vector<py::array_t<double>> quadratic_terms, const std::vector<py::array_t<int>> left_edges, const std::vector<py::array_t<int>> right_edges) {
                auto assignment_masks = get_assignment_masks(l, assignments);
                auto quadratic_masks = get_quadratic_terms_masks(l, quadratic_terms, left_edges, right_edges);
                return std::make_pair(assignment_masks, quadratic_masks);
                })
        .def("result_masks", [](const LPMP::multigraph_matching_input::labeling& l, const std::vector<py::array_t<double>> assignments, const std::vector<py::array_t<double>> quadratic_terms, const std::vector<py::array_t<int>> edges) {
                auto assignment_masks = get_assignment_masks(l, assignments);
                auto quadratic_masks = get_quadratic_terms_mask(l, quadratic_terms, edges);
                return std::make_pair(assignment_masks, quadratic_masks);
                })
        .def("__len__", [](const LPMP::multigraph_matching_input::labeling& l) { return l.size(); });


    py::class_<LPMP::multigraph_matching_input>(m, "multigraph_matching_input")
        .def(py::init<>())
        .def(py::init([](const std::vector<py::array_t<double>> assignments, const std::vector<py::array_t<double>> quadratic_terms, const std::vector<py::array_t<int>> left_edges, const std::vector<py::array_t<int>> right_edges) {
                LPMP::multigraph_matching_input instance;
                construct_from_arrays(instance, assignments, quadratic_terms, left_edges, right_edges); 
                return instance;
                }))
        .def(py::init([](const std::vector<py::array_t<double>> assignments, const std::vector<py::array_t<double>> quadratic_terms, const std::vector<py::array_t<int>> edges) {
                LPMP::multigraph_matching_input instance;
                construct_from_arrays(instance, assignments, quadratic_terms, edges); 
                return instance;
                }))
        .def("write", [](const LPMP::multigraph_matching_input& i, const std::string& filename) {
                std::fstream f;
                f.open(filename, std::ios::out);
                if(!f.is_open())
                    throw std::runtime_error("file " + filename + " could not be opened for exporting multigraph matching problem");
                return i.write_torresani_et_al(f); })
        .def("write", [](const LPMP::multigraph_matching_input& i){ return i.write_torresani_et_al(std::cout); })
        .def("evaluate", &LPMP::multigraph_matching_input::evaluate)
        ;


    using mgm_solver = LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::FMC_MGM<>>,LPMP::StandardTighteningVisitor>>; 
    py::class_<mgm_solver>(m, "multigraph_matching_message_passing_solver")
        .def(py::init<std::vector<std::string>&>())
        .def("construct", [](mgm_solver& s, const LPMP::multigraph_matching_input& input){ return s.GetProblemConstructor().construct(input); })
        .def("solve", &mgm_solver::Solve)
        .def("export_multigraph_matching_input", [](mgm_solver& s){ return s.GetProblemConstructor().export_multigraph_matching_input(); })
        .def("duality_gap", [](mgm_solver& s){
                const double lb = s.lower_bound();
                const double p = s.primal_cost();
                return p - lb;
                })
        .def("result",  [](mgm_solver& s){ return s.GetProblemConstructor().write_out_labeling(); });


}


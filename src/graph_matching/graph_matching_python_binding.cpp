#include "graph_matching/matching_problem_input.h"
#include "graph_matching/graph_matching.h"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"
#include <vector>
#include <tuple>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

void add_assignments(LPMP::graph_matching_input& instance, const py::array_t<double> assignments)
{
    const auto buf = assignments.request();
    const double* ptr = (double *) buf.ptr;
    const std::size_t dim1 = buf.shape[0];
    const std::size_t dim2 = buf.shape[1];
    for(std::size_t i=0; i<dim1; ++i) {
        for(std::size_t j=0; j<dim1; ++j) {
            //std::cout << "add assignment " << i*dim2+j << " = (" << i << "," << j << "), cost = " << ptr[i*dim2 + j] << "\n";
            instance.add_assignment(i, j, ptr[i*dim2 + j]);
        }
    }
}

void add_quadratic_terms(LPMP::graph_matching_input& instance, const py::array_t<double> quadratic_terms, const py::array_t<int> left_edges, const py::array_t<int> right_edges)
{
    const auto quadratic_buf = quadratic_terms.request();
    const auto left_edges_buf = left_edges.request();
    const auto right_edges_buf = right_edges.request();
    const double* quadratic_ptr = (double *) quadratic_buf.ptr;
    const int* left_edges_ptr = (int *) left_edges_buf.ptr;
    const int* right_edges_ptr = (int *) right_edges_buf.ptr;
    const std::size_t no_left_edges = quadratic_buf.shape[0];
    const std::size_t no_right_edges = quadratic_buf.shape[1];
    //std::cout << left_edges_buf.shape[0] << ", " << left_edges_buf.shape[1]  << "\n";
    assert(left_edges_buf.shape[1] == 2);
    assert(left_edges_buf.shape[0] == no_left_edges);
    //std::cout << right_edges_buf.shape[0] << ", " << right_edges_buf.shape[1]  << "\n";
    assert(right_edges_buf.shape[1] == 2);
    assert(right_edges_buf.shape[0] == no_right_edges);

    //std::cout << "nr left edges = " << no_left_edges << "\n";
    //std::cout << "nr right edges = " << no_right_edges << "\n";

    for(std::size_t i=0; i<no_left_edges; ++i) {
        for(std::size_t j=0; j<no_right_edges; ++j) {
            const std::size_t left_node_1 = left_edges_ptr[i*2];
            const std::size_t right_node_1 = right_edges_ptr[j*2];
            const std::size_t assignment_1 = left_node_1*instance.no_right_nodes + right_node_1;
            const std::size_t left_node_2 = left_edges_ptr[i*2+1];
            const std::size_t right_node_2 = right_edges_ptr[j*2+1];
            const std::size_t assignment_2 = left_node_2*instance.no_right_nodes + right_node_2;
            //std::cout << "left_edges[" << i << "] = (" << left_node_1 << "," << left_node_2 << ")\n";
            //std::cout << "assignment 1 = " << assignment_1 << "\n";
            //std::cout << "right_edges[" << j << "] = (" << right_node_1 << "," << right_node_2 << ")\n";
            //std::cout << "assignment 2 = " << assignment_2 << "\n";
            //std::cout << "quadratic term = " << quadratic_ptr[i*no_right_edges + j] << "\n";

            assert(instance.assignments[assignment_1].left_node == left_node_1);
            assert(instance.assignments[assignment_1].right_node == right_node_1);

            assert(instance.assignments[assignment_2].left_node == left_node_2);
            assert(instance.assignments[assignment_2].right_node == right_node_2);

            assert(left_node_1 != left_node_2);
            assert(right_node_1 != right_node_2);

            instance.add_quadratic_term(assignment_1, assignment_2, quadratic_ptr[i*no_right_edges + j]);
        }
    }
}

void construct_from_arrays(LPMP::graph_matching_input& instance, const py::array_t<double> assignments, const py::array_t<double> quadratic_terms, const py::array_t<int> left_edges, const py::array_t<int> right_edges)
{
    add_assignments(instance, assignments);
    add_quadratic_terms(instance, quadratic_terms, left_edges, right_edges);
}

py::array_t<char> get_assignment_mask(LPMP::graph_matching_input::labeling labeling, const py::array_t<double> assignments)
{
    char* assignment_mask = new char[assignments.size()];
    std::fill(assignment_mask, assignment_mask + assignments.size(), 0);
    assert(labeling.size() == assignments.shape(0));
    for(std::size_t i=0; i<labeling.size(); ++i) {
        if(labeling[i] != LPMP::linear_assignment_problem_input::no_assignment) {
            assert(labeling[i] < assignments.shape(1));
            assignment_mask[i*assignments.shape(0) + labeling[i]] = 1;
        }
    }
    return py::array({assignments.shape(0), assignments.shape(1)}, assignment_mask);
}

py::array_t<char> get_quadratic_terms_mask(LPMP::graph_matching_input::labeling labeling, const py::array_t<double> quadratic_terms, const py::array_t<int> left_edges, const py::array_t<int> right_edges)
{
    const auto l = left_edges.unchecked<2>();
    const auto r = right_edges.unchecked<2>();
    char* quadratic_mask = new char[quadratic_terms.size()];
    std::fill(quadratic_mask, quadratic_mask + quadratic_terms.size(), 0);
    for(std::size_t i=0; i<quadratic_terms.shape(0); ++i) {
        for(std::size_t j=0; j<quadratic_terms.shape(1); ++j) {
            const std::size_t left_node_1 = l(i,0);
            const std::size_t right_node_1 = r(j,0);
            const std::size_t left_node_2 = l(i,1);
            const std::size_t right_node_2 = r(j,1);

            if(labeling[left_node_1] == right_node_1 && labeling[left_node_2] == right_node_2)
                quadratic_mask[i*quadratic_terms.shape(0) + j] = 1;
        }
    } 

    return py::array({quadratic_terms.shape(0), quadratic_terms.shape(1)}, quadratic_mask);
}



PYBIND11_MODULE(graph_matching_py, m) {
    m.doc() = "python binding for LPMP graph matching";

    py::class_<LPMP::graph_matching_input::labeling>(m, "graph_matching_labeling")
        .def(py::init<>())
        .def(py::init([](const std::vector<std::size_t>& labels) {
                    LPMP::graph_matching_input::labeling l;
                    for(const auto x : labels) l.push_back(x);
                    return l;
                    }))
        .def("size", &LPMP::graph_matching_input::labeling::size)
        .def("__getitem__", [](const LPMP::graph_matching_input::labeling &l, size_t i) { 
                if (i >= l.size()) throw py::index_error();
                return l[i];
                })
        .def("result_masks", [](const LPMP::graph_matching_input::labeling& l, const py::array_t<double> assignments, const py::array_t<double> quadratic_terms, const py::array_t<std::size_t> left_edges, const py::array_t<std::size_t> right_edges) {
                auto assignment_mask = get_assignment_mask(l, assignments);
                auto quadratic_mask = get_quadratic_terms_mask(l, quadratic_terms, left_edges, right_edges);
                return std::make_pair(assignment_mask, quadratic_mask);
                });

    py::class_<LPMP::linear_assignment_problem_input>(m, "linear_assignment_problem_input")
        .def(py::init<>())
        .def("add_assignment", &LPMP::linear_assignment_problem_input::add_assignment)
        .def("evaluate", &LPMP::linear_assignment_problem_input::evaluate);

    m.def("graph_matching_no_assignment", []() { return LPMP::linear_assignment_problem_input::no_assignment; });


    py::class_<LPMP::graph_matching_input, LPMP::linear_assignment_problem_input>(m, "graph_matching_input")
        .def(py::init<>())
        .def(py::init([](const py::array_t<double> assignments, const py::array_t<double> quadratic_terms, const py::array_t<std::size_t> left_edges, const py::array_t<std::size_t> right_edges) {
                LPMP::graph_matching_input instance;
                construct_from_arrays(instance, assignments, quadratic_terms, left_edges, right_edges); 
                return instance;
                }))
        .def("add_assignment", &LPMP::graph_matching_input::add_assignment)
        //.def("add_assignment", &LPMP::linear_assignment_problem_input::add_assignment)
        .def("add_quadratic_term", &LPMP::graph_matching_input::add_quadratic_term)
        .def("write", [](const LPMP::graph_matching_input& i){ return i.write(std::cout); })
        .def("evaluate", &LPMP::graph_matching_input::evaluate)
        .def("add_assignments", add_assignments)
        .def("add_quadratic_terms", add_quadratic_terms);


    using gm_mp_solver = LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::FMC_MP>,LPMP::StandardVisitor>>; 
    py::class_<gm_mp_solver>(m, "graph_matching_message_passing_solver")
        .def(py::init<std::vector<std::string>&>())
        .def("construct", [](gm_mp_solver& s, const LPMP::graph_matching_input& input){ return s.GetProblemConstructor().construct(input); })
        .def("solve", &gm_mp_solver::Solve)
        .def("export", [](gm_mp_solver& s){ return s.GetProblemConstructor().export_graph_matching_input(); })
        .def("result", [](gm_mp_solver& s){ return s.GetProblemConstructor().write_out_labeling(); });

    using gm_mp_q_solver = LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::FMC_MP_Q>,LPMP::StandardVisitor>>; 
    py::class_<gm_mp_q_solver>(m, "graph_matching_message_passing_interquadratic_message_solver")
        .def(py::init<std::vector<std::string>&>())
        .def("construct", [](gm_mp_q_solver& s, const LPMP::graph_matching_input& input){ return s.GetProblemConstructor().construct(input); })
        .def("solve", &gm_mp_q_solver::Solve)
        .def("export",  [](gm_mp_q_solver& s){ return s.GetProblemConstructor().export_graph_matching_input(); })
        .def("result",  [](gm_mp_q_solver& s){ return s.GetProblemConstructor().write_out_labeling(); });

    using gm_mp_t_solver = LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::FMC_MP_T>,LPMP::StandardVisitor>>; 
    py::class_<gm_mp_t_solver>(m, "graph_matching_message_passing_tightening_solver")
        .def(py::init<std::vector<std::string>&>())
        .def("construct", [](gm_mp_t_solver& s, const LPMP::graph_matching_input& input){ return s.GetProblemConstructor().construct(input); })
        .def("solve", &gm_mp_t_solver::Solve)
        .def("export",  [](gm_mp_t_solver& s){ return s.GetProblemConstructor().export_graph_matching_input(); })
        .def("result",  [](gm_mp_t_solver& s){ return s.GetProblemConstructor().write_out_labeling(); });

    using gm_mp_q_t_solver = LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::FMC_MP_Q_T>,LPMP::StandardVisitor>>; 
    py::class_<gm_mp_q_t_solver>(m, "graph_matching_message_passing_interquadratic_message_tightening_solver")
        .def(py::init<std::vector<std::string>&>())
        .def("construct", [](gm_mp_q_t_solver& s, const LPMP::graph_matching_input& input){ return s.GetProblemConstructor().construct(input); })
        .def("solve", &gm_mp_q_t_solver::Solve)
        .def("export",  [](gm_mp_q_t_solver& s){ return s.GetProblemConstructor().export_graph_matching_input(); })
        .def("result",  [](gm_mp_q_t_solver& s){ return s.GetProblemConstructor().write_out_labeling(); });

    using gm_mrf_solver = LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::FMC_GM>,LPMP::StandardVisitor>>; 
    py::class_<gm_mrf_solver>(m, "graph_matching_mrf_solver")
        .def(py::init<std::vector<std::string>&>())
        .def("construct", [](gm_mrf_solver& s, const LPMP::graph_matching_input& input){ return s.GetProblemConstructor().construct(input); })
        .def("solve", &gm_mrf_solver::Solve)
        .def("export",  [](gm_mrf_solver& s){ return s.GetProblemConstructor().export_graph_matching_input(); })
        .def("result",  [](gm_mrf_solver& s){ return s.GetProblemConstructor().write_out_labeling(); });

    using gm_mrf_t_solver = LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::FMC_GM_T>,LPMP::StandardVisitor>>; 
    py::class_<gm_mrf_t_solver>(m, "graph_matching_mrf_tightening_solver")
        .def(py::init<std::vector<std::string>&>())
        .def("construct", [](gm_mrf_t_solver& s, const LPMP::graph_matching_input& input){ return s.GetProblemConstructor().construct(input); })
        .def("solve", &gm_mrf_t_solver::Solve)
        .def("export",  [](gm_mrf_t_solver& s){ return s.GetProblemConstructor().export_graph_matching_input(); })
        .def("result",  [](gm_mrf_t_solver& s){ return s.GetProblemConstructor().write_out_labeling(); });

}

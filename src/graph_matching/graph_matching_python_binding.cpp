#include "graph_matching/matching_problem_input.h"
#include "graph_matching/graph_matching.h"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(graph_matching_python_binding, m) {
    m.doc() = "python binding for LPMP graph matching";

    py::class_<LPMP::graph_matching_input::labeling>(m, "graph_matching_labeling")
        .def(py::init<>());

    py::class_<LPMP::graph_matching_input>(m, "graph_matching_input")
        .def(py::init<>())
        .def("add_assignment", &LPMP::graph_matching_input::add_assignment)
        .def("add_quadratic_term", &LPMP::graph_matching_input::add_quadratic_term)
        .def("evaluate", &LPMP::graph_matching_input::evaluate);

    using gm_mp_solver = LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::FMC_MP>,LPMP::StandardVisitor>>; 
    py::class_<gm_mp_solver>(m, "graph_matching_message_passing_solver")
        .def(py::init<std::vector<std::string>&>())
        .def("construct", [](gm_mp_solver& s, const LPMP::graph_matching_input& input){ return s.GetProblemConstructor().construct(input); })
        .def("solve", &gm_mp_solver::Solve)
        .def("result",  [](gm_mp_solver& s){ return s.GetProblemConstructor().write_out_labeling(); });

    using gm_mp_q_solver = LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::FMC_MP_Q>,LPMP::StandardVisitor>>; 
    py::class_<gm_mp_q_solver>(m, "graph_matching_message_passing_interquadratic_message_solver")
        .def(py::init<std::vector<std::string>&>())
        .def("construct", [](gm_mp_q_solver& s, const LPMP::graph_matching_input& input){ return s.GetProblemConstructor().construct(input); })
        .def("solve", &gm_mp_q_solver::Solve)
        .def("result",  [](gm_mp_q_solver& s){ return s.GetProblemConstructor().write_out_labeling(); });

    using gm_mp_t_solver = LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::FMC_MP_T>,LPMP::StandardVisitor>>; 
    py::class_<gm_mp_t_solver>(m, "graph_matching_message_passing_tightening_solver")
        .def(py::init<std::vector<std::string>&>())
        .def("construct", [](gm_mp_t_solver& s, const LPMP::graph_matching_input& input){ return s.GetProblemConstructor().construct(input); })
        .def("solve", &gm_mp_t_solver::Solve)
        .def("result",  [](gm_mp_t_solver& s){ return s.GetProblemConstructor().write_out_labeling(); });

    using gm_mp_q_t_solver = LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::FMC_MP_Q_T>,LPMP::StandardVisitor>>; 
    py::class_<gm_mp_q_t_solver>(m, "graph_matching_message_passing_interquadratic_message_tightening_solver")
        .def(py::init<std::vector<std::string>&>())
        .def("construct", [](gm_mp_q_t_solver& s, const LPMP::graph_matching_input& input){ return s.GetProblemConstructor().construct(input); })
        .def("solve", &gm_mp_q_t_solver::Solve)
        .def("result",  [](gm_mp_q_t_solver& s){ return s.GetProblemConstructor().write_out_labeling(); });

    using gm_mrf_solver = LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::FMC_GM>,LPMP::StandardVisitor>>; 
    py::class_<gm_mrf_solver>(m, "graph_matching_mrf_solver")
        .def(py::init<std::vector<std::string>&>())
        .def("construct", [](gm_mrf_solver& s, const LPMP::graph_matching_input& input){ return s.GetProblemConstructor().construct(input); })
        .def("solve", &gm_mrf_solver::Solve)
        .def("result",  [](gm_mrf_solver& s){ return s.GetProblemConstructor().write_out_labeling(); });

    using gm_mrf_t_solver = LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::FMC_GM_T>,LPMP::StandardVisitor>>; 
    py::class_<gm_mrf_t_solver>(m, "graph_matching_mrf_tightening_solver")
        .def(py::init<std::vector<std::string>&>())
        .def("construct", [](gm_mrf_t_solver& s, const LPMP::graph_matching_input& input){ return s.GetProblemConstructor().construct(input); })
        .def("solve", &gm_mrf_t_solver::Solve)
        .def("result",  [](gm_mrf_t_solver& s){ return s.GetProblemConstructor().write_out_labeling(); });

}

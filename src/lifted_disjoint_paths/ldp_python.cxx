#include <pybind11/pybind11.h>
//#include "disjoint-paths/disjointPathsMethods.hxx"
#include "lifted_disjoint_paths/ldp_config.hxx"
#include "visitors/standard_visitor.hxx"
#include "lifted_disjoint_paths/ldp_instance.hxx"
#include "solver.hxx"
#include "lifted_disjoint_paths/lifted_disjoint_paths_fmc.h"
#include "LP.h"




namespace py = pybind11;



PYBIND11_MODULE(ldpMessagePassing, m) {
     m.doc() = "python binding for lifted disjoint paths based on message passing";

     py::class_<disjointPaths::ParametersParser>(m,"ParametersParser")
             .def(py::init<>())
             .def("get_parsed_params",&disjointPaths::ParametersParser::getParsedStrings,"getting the parsed strings from parser")
             .def("init_from_file", py::overload_cast<std::string&>(&disjointPaths::ParametersParser::initFromFile),"Parses parameters from a file");

     py::class_<LPMP::lifted_disjoint_paths::ConfigDisjoint<size_t>>(m, "LdpParams")
        .def(py::init<std::map<std::string,std::string>&>());

     py::class_<disjointPaths::VertexGroups<>>(m, "TimeFramesToVertices")
        .def(py::init<>())
        .def("init_from_vector", &disjointPaths::VertexGroups<>::initFromVector, "Initializes vertices in time frames from a vector of size_t")
        .def("init_from_file", &disjointPaths::VertexGroups<>::initFromFile<LPMP::lifted_disjoint_paths::ConfigDisjoint<>>, "Initializes vertices in time frames from a file");

     py::class_<disjointPaths::CompleteStructure<>>(m, "GraphStructure")
        .def(py::init<disjointPaths::VertexGroups<> &>())
        .def("add_edges_from_array", &disjointPaths::CompleteStructure<>::addEdgesFromMatrix, "Initializes edges of the graph between two time frames from a matrix.")
        .def("add_edges_from_file", &disjointPaths::CompleteStructure<>::addEdgesFromFile<LPMP::lifted_disjoint_paths::ConfigDisjoint<>>, "Initializes all edges of the graph from a file.");

     py::class_<LPMP::lifted_disjoint_paths::LdpInstance>(m, "LdpInstance")
        .def(py::init<LPMP::lifted_disjoint_paths::ConfigDisjoint<size_t> &,disjointPaths::CompleteStructure<>&>());

     py::class_<LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::lifted_disjoint_paths_FMC>,LPMP::StandardTighteningVisitor>>>(m,"Solver")
             .def(py::init<std::vector<std::string>&>())
             .def("solve",&LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::lifted_disjoint_paths_FMC>,LPMP::StandardTighteningVisitor>>::Solve);


     m.def("construct",&LPMP::constructProblemFromSolver<LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::lifted_disjoint_paths_FMC>,LPMP::StandardTighteningVisitor>>,LPMP::lifted_disjoint_paths::LdpInstance>,"constructing problem from instance");




    // m.def("solve_ilp", py::overload_cast<disjointPaths::DisjointParams<>&, disjointPaths::CompleteStructure<>&>(&disjointPaths::solver_ilp_intervals<>), "Solve lifted disjoint paths");

   //  m.def("solve_message_passing", &disjointPaths::solver_ilp<>, "Solve lifted disjoint paths");

     //m.def("write_output_to_file", &disjointPaths::writeOutputToFile, "Write output tracks to a specified file");





}


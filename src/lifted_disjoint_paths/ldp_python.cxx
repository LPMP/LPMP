#include <pybind11/pybind11.h>
//#include "disjoint-paths/disjointPathsMethods.hxx"
#include "lifted_disjoint_paths/ldp_parameters.hxx"
#include "visitors/standard_visitor.hxx"
#include "lifted_disjoint_paths/ldp_instance.hxx"
#include "solver.hxx"
#include "lifted_disjoint_paths/lifted_disjoint_paths_fmc.h"
#include "LP.h"
#include "lifted_disjoint_paths/ldp_interval_connection.hxx"
#include "lifted_disjoint_paths/ldp_paths_extractor.hxx"



namespace py = pybind11;



PYBIND11_MODULE(ldpMessagePassingPy, m) {
     m.doc() = "python binding for lifted disjoint paths based on message passing";

     py::class_<LPMP::ParametersParser>(m,"ParametersParser")
             .def(py::init<>())
             .def("get_parsed_params",&LPMP::ParametersParser::getParsedStrings,"getting the parsed strings from parser")
             .def("init_from_file", py::overload_cast<std::string&>(&LPMP::ParametersParser::initFromFile),"Parses parameters from a file");

     py::class_<LPMP::lifted_disjoint_paths::LdpParameters<size_t>>(m, "LdpParams")
        .def(py::init<std::map<std::string,std::string>&>());

     py::class_<LPMP::VertexGroups<>>(m, "TimeFramesToVertices")
        .def(py::init<>())
        .def("get_vertex_shift",&LPMP::VertexGroups<>::getVertexShiftBack,"In case that the object was initialized from a file with minimal time frame higher than one, this gives the ID of the minimal vertex in the minimal used time layer.")
        .def("init_from_vector", py::overload_cast<const std::vector<size_t>&>(&LPMP::VertexGroups<>::initFromVector), "Initializes vertices in time frames from a vector of size_t")
        .def("init_from_vector_shift", py::overload_cast<const std::vector<size_t>&,size_t ,size_t>(&LPMP::VertexGroups<>::initFromVector), "Initializes vertices in time frames from a vector of size_t,time shift and vertex shift")
        .def("get_max_time",&LPMP::VertexGroups<>::getMaxTime,"returns maximal time")
        .def("init_from_file", &LPMP::VertexGroups<>::initFromFile<LPMP::lifted_disjoint_paths::LdpParameters<>>, "Initializes vertices in time frames from a file");


     py::class_<LPMP::CompleteStructure<>>(m, "GraphStructure")
             .def(py::init<LPMP::VertexGroups<> &>())
            // .def("add_edges_from_array", &LPMP::CompleteStructure<>::addEdgesFromMatrix, "Initializes edges of the graph between two time frames from a matrix.")
             .def("add_edges_from_file", &LPMP::CompleteStructure<>::addEdgesFromFile<LPMP::lifted_disjoint_paths::LdpParameters<>>, "Initializes all edges of the graph from a file.")
             .def("add_edges_from_vectors", &LPMP::CompleteStructure<>::addEdgesFromVectors<LPMP::lifted_disjoint_paths::LdpParameters<>>, "Initializes edges of the graph from an Nx2 array of size_t with edge vertices and an Nx1 array of doubles with costs. Restrictions on maximal vertex and maximal time gap from parameters apply.")
             .def("add_edges_from_vectors_all", &LPMP::CompleteStructure<>::addEdgesFromVectorsAll, "Initializes edges of the graph from an Nx2 array of size_t with edge vertices and an Nx1 array of doubles with costs. All edges added. No restriction on maximal timegap. ")
             //.def("get_edge_labels",&LPMP::CompleteStructure<>::getGraphEdgeLabels,"Returns 0/1 labels of all input edges w.r.t. given set of paths. Label one is given iff detections belong to the same path." )
             .def("set_score_of_vertices",&LPMP::CompleteStructure<>::setVerticesCosts,"Expects array with score of all graph vertices.")
             .def("get_edge_list",&LPMP::CompleteStructure<>::getEdgeList,"Return list of edges present in this graph structure.");

     py::class_<LPMP::lifted_disjoint_paths::LdpInstance>(m, "LdpInstance")
             .def(py::init<LPMP::lifted_disjoint_paths::LdpParameters<size_t> &,LPMP::CompleteStructure<>&>())
             .def(py::init<LPMP::lifted_disjoint_paths::LdpParameters<size_t>&,LPMP::LdpBatchProcess&>())
             .def(py::init<LPMP::lifted_disjoint_paths::LdpParameters<size_t>&,LPMP::LdpIntervalConnection&>())
             .def(py::init<LPMP::lifted_disjoint_paths::LdpParameters<>&,const py::array_t<size_t>&,const py::array_t<size_t>&,const  py::array_t<double>& ,const  py::array_t<double>&,const  py::array_t<double>&,LPMP::VertexGroups<>&>()) ;


     using problemSolver=LPMP::ProblemConstructorRoundingSolver<LPMP::Solver<LPMP::LP<LPMP::lifted_disjoint_paths_FMC>,LPMP::StandardTighteningVisitor>>;
     py::class_<problemSolver>(m,"Solver")
             .def(py::init<std::vector<std::string>&>())
             .def("solve",&problemSolver::Solve)
             .def("get_lower_bound",&problemSolver::lower_bound,"Returns lower bound")
             .def("get_best_primal_value",&problemSolver::primal_cost,"returns best primal value")
             .def("get_best_primal", [](problemSolver &solver) {return solver.GetProblemConstructor().getBestPrimal(); },"Returns paths obtained from best so far primal solution.");





     m.def("construct",&LPMP::constructProblemFromSolver<problemSolver,LPMP::lifted_disjoint_paths::LdpInstance>,"constructing problem from instance");


     m.def("get_base_edge_labels",&LPMP::getBaseEdgeLabels<std::vector<std::array<size_t,2>>>,"Given a vector of base edge vertices, vector of solution paths and the number of graph vertices, it returns labels to base edges.");

     m.def("get_lifted_edge_labels",&LPMP::getLiftedEdgeLabels<std::vector<std::array<size_t,2>>>,"Given a vector of lifted edge vertices, vector of solution paths and the number of graph vertices, it returns labels to base edges.");


     py::class_<LPMP::LdpBatchProcess>(m,"BatchProcess")
             .def(py::init<LPMP::VertexGroups<>&, std::vector<std::array<size_t,2>>&, size_t, size_t, size_t, size_t>(),"TimeFramesToVertices,n x 2 vector of vertex ID->label, max used label, max time for using labels,min time of batch, max time of batch")
             .def("init_edges_from_vectors",&LPMP::LdpBatchProcess::initEdgesFromVector,"Requires n x 2 vector of edge vertices and n x 1 vector of costs")
             .def("init_vertices_from_vectors",&LPMP::LdpBatchProcess::initVertexScoreFromVector,"Requires n x 1 list of vertices and n x 1 list of their costs")
             .def("init_from_file",&LPMP::LdpBatchProcess::initFromFile,"Initializes both edges and vectors from a file")
             .def("decode_solution",&LPMP::LdpBatchProcess::decode,"Given paths resulting from solver, creates labels for new vertices in form: vector n x 2: vertexID->label")
             .def("get_labels",&LPMP::LdpBatchProcess::getDecodedLabels,"Returns labels obtained from solution in form: vector n x 2: vertexID->label")
             .def("get_index_to_delete",&LPMP::LdpBatchProcess::getIndexToDel,"Returns index than needs to be used for deleting outdated labels in label vector")
             .def("get_max_used_label",&LPMP::LdpBatchProcess::getMaxLabelsSoFar,"Returns max used label, needed for constructor of next batch");

     py::class_<LPMP::LdpPathsExtractor>(m,"PathsExtractor")
             .def(py::init<LPMP::VertexGroups<>&,const std::vector<std::vector<size_t>>&, size_t,bool,bool,size_t>(),"Requires: TimeFramesToVertices of the respective interval, solution paths,length to cut off, bool is this first interval,bool is this last interval,vertex shift (global id of the first valid vertex)")
             .def("get_extracted_paths",&LPMP::LdpPathsExtractor::getExtractedPaths,"Returns list of paths. Where vertices have their global IDs");

//     py::class_<LPMP::LdpIntervalConnectionResult>(m,"IntervalConnectionResults")
//             .def(py::init<>())
//             .def("get_middle_paths",&LPMP::LdpIntervalConnectionResult::getMiddleToSecond,"Getting paths between intervals.")
//             .def("get_start_in_middle",&LPMP::LdpIntervalConnectionResult::getStartInMiddle,"Getting paths in the middle interval that starting.")
//             .def("get_start_in_second",&LPMP::LdpIntervalConnectionResult::getStartInSecond,"Getting paths from the second interval that are starting")
//             .def("get_middle_to_second",&LPMP::LdpIntervalConnectionResult::getMiddleToSecond,"Mapping of middle interval paths to second interval paths.")
//             .def("get_first_to_middle",&LPMP::LdpIntervalConnectionResult::getFirstToMiddle,"Mapping of first interval paths to middle interval paths.");

     py::class_<LPMP::LdpLabelsAssignment>(m,"LabelAssignment")
             .def(py::init<>())
             .def("init",&LPMP::LdpLabelsAssignment::init,"Initializes label assignment with paths from the first interval.")
             .def("get_labels",&LPMP::LdpLabelsAssignment::getVerticesToLabels,"Returns vector nx2 of pairs: vertex ID-> vertex label")
             .def("update",&LPMP::LdpLabelsAssignment::update,"Updates labels given new paths and connections from old to new paths.");

     py::class_<LPMP::LdpIntervalConnection>(m,"IntervalConnection")
             .def(py::init<const LPMP::LdpPathsExtractor&, const LPMP::LdpPathsExtractor&>(),"Requires: PathsExtractors of the first and of the second interval")
             .def("decode_paths",&LPMP::LdpIntervalConnection::decodePaths,"Decodes paths returned from the solver into the global vertex IDs.")
             .def("init_from_file",&LPMP::LdpIntervalConnection::initFromFile,"Initializes both edges and vectors from a file")
             .def("init_vertices_from_vectors",&LPMP::LdpIntervalConnection::initScoreOfVertices,"Requires n x 1 list of vertices and n x 1 list of their costs. Vertices not provided in the list have cost zero.")
             .def("init_edges_from_vectors",&LPMP::LdpIntervalConnection::initEdgesFromVectors,"Requires n x 2 vector of edge vertices and n x 1 vector of costs for edges in the two intervals. Extra edges are ignored.")
             .def("get_middle_paths",&LPMP::LdpIntervalConnection::getMiddleIntervalPaths,"Getting paths between intervals.")
             .def("get_middle_to_second",&LPMP::LdpIntervalConnection::getMiddleToSecond,"Mapping of middle interval paths to second interval paths.")
             .def("create_result_structures",&LPMP::LdpIntervalConnection::createResultsStructures,"Decodes output paths and creates resulting structures.")
             .def("get_first_to_middle",&LPMP::LdpIntervalConnection::getFirstToMiddle,"Mapping of first interval paths to middle interval paths.");






}


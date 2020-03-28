#include "graph_matching/graph_matching_python_binding_helper.h"

namespace LPMP {
    namespace py_helper {

    void add_assignments(LPMP::graph_matching_input& instance, const py::array_t<double> assignments)
    {
        const auto _assignments = assignments.unchecked<2>();
        const std::size_t dim1 = _assignments.shape(0);
        const std::size_t dim2 = _assignments.shape(1);
        for(std::size_t i=0; i<dim1; ++i) {
            for(std::size_t j=0; j<dim2; ++j) {
                //std::cout << "add assignment " << i*dim2+j << " = (" << i << "," << j << "), cost = " << _assignments(i,j) << "\n";
                instance.add_assignment(i, j, _assignments(i,j));
            }
        }
    }

    void add_quadratic_terms(LPMP::graph_matching_input& instance, const py::array_t<double> quadratic_terms, const py::array_t<int> left_edges, const py::array_t<int> right_edges)
    {
        const auto _left_edges = left_edges.unchecked<2>();
        const auto _right_edges = right_edges.unchecked<2>();
        const auto _quadratic_terms = quadratic_terms.unchecked<2>();

        const std::size_t no_left_edges = quadratic_terms.shape(0);
        const std::size_t no_right_edges = quadratic_terms.shape(1);
        //std::cout << no_left_edges << ";" << left_edges.shape(0) << ", " << left_edges.shape(1)  << "\n";
        if(left_edges.shape(1) != 2)
            throw std::runtime_error("graph matching python binding: left_edges must have second dimension = 2");
        if(left_edges.shape(0) != no_left_edges)
            throw std::runtime_error("graph matching python binding: dimension incompatibility");
        //std::cout << right_edges.shape(0) << ", " << right_edges.shape(1)  << "\n";
        if(right_edges.shape(1) != 2)
            throw std::runtime_error("graph matching python binding: right must have second dimension = 2");
        if(right_edges.shape(0) != no_right_edges)
            throw std::runtime_error("graph matching python binding: dimension incompatibility");

        //std::cout << "nr left edges = " << no_left_edges << "\n";
        //std::cout << "nr right edges = " << no_right_edges << "\n";

        for(std::size_t i=0; i<no_left_edges; ++i) {
            for(std::size_t j=0; j<no_right_edges; ++j) {
                const std::size_t left_node_1 = _left_edges(i,0);//left_edges_ptr[i*2];
                const std::size_t right_node_1 = _right_edges(j,0);//right_edges_ptr[j*2];
                const std::size_t assignment_1 = left_node_1*instance.no_right_nodes + right_node_1;
                if(left_node_1 >= instance.no_left_nodes)
                    throw std::runtime_error("left node out of bounds");
                if(right_node_1 >= instance.no_right_nodes) {
                    throw std::runtime_error("right node out of bounds");
                }

                const std::size_t left_node_2 = _left_edges(i,1);//left_edges_ptr[i*2+1];
                const std::size_t right_node_2 = _right_edges(j,1);//right_edges_ptr[j*2+1];
                const std::size_t assignment_2 = left_node_2*instance.no_right_nodes + right_node_2;
                if(left_node_2 >= instance.no_left_nodes)
                    throw std::runtime_error("left node out of bounds");
                if(right_node_2 >= instance.no_right_nodes) {
                    throw std::runtime_error("right node out of bounds");
                }

                //std::cout << "left_edges[" << i << "] = (" << left_node_1 << "," << left_node_2 << ")\n";
                //std::cout << "assignment 1 = " << assignment_1 << "\n";
                //std::cout << "right_edges[" << j << "] = (" << right_node_1 << "," << right_node_2 << ")\n";
                //std::cout << "assignment 2 = " << assignment_2 << "\n";
                //std::cout << "quadratic term = " << _quadratic_terms(i,j) << "\n";

                assert(instance.assignments[assignment_1].left_node == left_node_1);
                assert(instance.assignments[assignment_1].right_node == right_node_1);

                assert(instance.assignments[assignment_2].left_node == left_node_2);
                assert(instance.assignments[assignment_2].right_node == right_node_2);

                if(left_node_1 == left_node_2) 
                    throw std::runtime_error("left nodes must be distinct");
                if(right_node_1 == right_node_2)
                    throw std::runtime_error("right nodes must be distinct");

                instance.add_quadratic_term(assignment_1, assignment_2, _quadratic_terms(i,j));
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
        assert(assignments.size() == assignments.shape(0)*assignments.shape(1));
        std::fill(assignment_mask, assignment_mask + assignments.size(), 0);
        if(labeling.size() != assignments.shape(0))
            throw std::runtime_error("labeling must be of equal size as first dimension of assignment matrix");
        for(std::size_t i=0; i<labeling.size(); ++i) {
            if(labeling[i] != LPMP::linear_assignment_problem_input::no_assignment) {
                if(labeling[i] >= assignments.shape(1))
                    throw std::runtime_error("labeling entry not valid");

                assert(i*assignments.shape(1) + labeling[i] < assignments.size());
                assignment_mask[i*assignments.shape(1) + labeling[i]] = 1;
            }
        }
        return py::array({assignments.shape(0), assignments.shape(1)}, assignment_mask);
    }

    py::array_t<char> get_quadratic_terms_mask(LPMP::graph_matching_input::labeling labeling, const py::array_t<double> quadratic_terms, const py::array_t<int> left_edges, const py::array_t<int> right_edges)
    {
        const auto l = left_edges.unchecked<2>();
        const auto r = right_edges.unchecked<2>();
        char* quadratic_mask = new char[quadratic_terms.size()];
        assert(quadratic_terms.size() == quadratic_terms.shape(0)*quadratic_terms.shape(1));
        std::fill(quadratic_mask, quadratic_mask + quadratic_terms.size(), 0);
        for(std::size_t i=0; i<quadratic_terms.shape(0); ++i) {
            for(std::size_t j=0; j<quadratic_terms.shape(1); ++j) {
                const std::size_t left_node_1 = l(i,0);
                const std::size_t right_node_1 = r(j,0);
                const std::size_t left_node_2 = l(i,1);
                const std::size_t right_node_2 = r(j,1);

                if(labeling[left_node_1] == right_node_1 && labeling[left_node_2] == right_node_2) {
                    assert(i*quadratic_terms.shape(1) + j < quadratic_terms.size());
                    quadratic_mask[i*quadratic_terms.shape(1) + j] = 1;
                }
            }
        } 

        return py::array({quadratic_terms.shape(0), quadratic_terms.shape(1)}, quadratic_mask);
    }

    }
}

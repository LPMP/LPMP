#include "graph_matching/graph_matching_instance.h"
#include <algorithm>
#include <array>

namespace LPMP {

    // LAP_instance

    void LAP_instance::add_assignment(const std::size_t left_node, const std::size_t right_node, const double cost)
    {
        assignments_.push_back({left_node, right_node, cost}); 
        no_left_nodes_ = std::max(no_left_nodes_, left_node+1);
        no_right_nodes_ = std::max(no_right_nodes_, right_node+1);
    }

    std::size_t LAP_instance::no_mcf_nodes() const { return no_left_nodes() + no_right_nodes() + 2; } 
    std::size_t LAP_instance::no_mcf_edges() const { return assignments().size() + no_left_nodes() + no_right_nodes() + 1; }
    std::size_t LAP_instance::no_left_nodes() const { return no_left_nodes_; }
    std::size_t LAP_instance::no_right_nodes() const { return no_right_nodes_; }

    void LAP_instance::normalize()
    {
        // merge parallel edges
        std::sort(assignments_.begin(), assignments_.end(), [](const auto& a1, const auto& a2) { 
                std::array<std::size_t,2> e1 {a1.left_node, a1.right_node};
                std::array<std::size_t,2> e2 {a2.left_node, a2.right_node};
                return std::lexicographical_compare(e1.begin(), e1.end(), e2.begin(), e2.end());
                });

        std::vector<assignment> normalized_assignments; 

        // merge matching edge copies and add them to edges
        for(std::size_t k=0; k<assignments_.size();) {
            normalized_assignments.push_back(assignments_[k]);
            ++k; 
            while(k<assignments_.size() && normalized_assignments.back().left_node == assignments_[k].left_node && normalized_assignments.back().right_node == assignments_[k].right_node) {
                normalized_assignments.back().cost += assignments_[k].cost;
                ++k; 
            }
        }

        std::swap(normalized_assignments, assignments_); 
    }

    double LAP_instance::evaluate(const LAP_labeling& l) const
    {
        if(l.no_left_nodes() != no_left_nodes())
            throw std::runtime_error("labeling must have equal number of left nodes as matching problem.");

        double cost = 0.0;

        std::vector<std::size_t> nodes_taken(no_right_nodes(),0);
        for(const auto& a : assignments()) {
            if(l[a.left_node] == a.right_node) {
                cost += a.cost;
                assert(a.left_node < nodes_taken.size());
                nodes_taken[a.left_node]++;
            }
        }

        if(*std::max_element(nodes_taken.begin(), nodes_taken.end()) > 1)
            return std::numeric_limits<double>::infinity();

        return cost;
    }

    // LAP_labeling

    bool LAP_labeling::check_primal_consistency() const
    {
        std::vector<char> labels_taken(this->highest_matched_node()+1, 0);
        for(const std::size_t l : *this) {
            if(l == std::numeric_limits<std::size_t>::max())
                continue;
            if(labels_taken[l] > 0)
                return false;
            labels_taken[l] = 1;
        }
        return true;
    }

    std::size_t LAP_labeling::highest_matched_node() const
    {
        std::size_t max_label = 0;
        for(std::size_t i=0; i<this->size(); ++i)
            if((*this)[i] != std::numeric_limits<std::size_t>::max())
                max_label = std::max(max_label, (*this)[i]);
        return max_label;
    }

    // graph_matching_instance

    double graph_matching_instance::evaluate(const LAP_labeling& l) const
    {
        double cost = LAP_instance::evaluate(l);
        for(const auto& q : quadratic_terms()) {
            const auto a1 = this->assignments()[q.assignment_1];
            const auto a2 = this->assignments()[q.assignment_2];
            if(l[a1.left_node] == a1.right_node && l[a2.left_node] == a2.right_node)
                cost += q.cost;
        }
        return cost;
    }
}

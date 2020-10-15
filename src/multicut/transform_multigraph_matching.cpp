#include <algorithm>
#include <numeric>
#include <cassert>
#include "multicut/transform_multigraph_matching.h"
#include "union_find.hxx"

namespace LPMP {

    multigraph_matching_correlation_clustering_transform::multigraph_matching_correlation_clustering_transform()
    {}

    multigraph_matching_correlation_clustering_transform::multigraph_matching_correlation_clustering_transform(std::shared_ptr<multigraph_matching_input> _mgm)
    {
        set_multigraph_matching_input(_mgm);
    }

    void multigraph_matching_correlation_clustering_transform::set_multigraph_matching_input(std::shared_ptr<multigraph_matching_input> _mgm)
    {
        cc = std::make_shared<correlation_clustering_instance>();
        mgm = _mgm;
        compute_no_nodes();
        compute_node_offsets();
        transform_multigraph_matching_to_correlation_clustering(); 
    }

    std::size_t multigraph_matching_correlation_clustering_transform::no_graphs() const 
    {
        return no_nodes_.size(); 
    }
    std::size_t multigraph_matching_correlation_clustering_transform::total_no_nodes() const 
    { 
        return graph_node_offsets_.back() + no_nodes_.back(); 
    }
    std::size_t multigraph_matching_correlation_clustering_transform::no_nodes(const std::size_t i) const 
    { 
        assert(i < no_graphs());
        return no_nodes_[i]; 
    }

    std::size_t multigraph_matching_correlation_clustering_transform::multigraph_matching_to_correlation_clustering_node(const std::size_t graph_no, const std::size_t node_no) const
    {
        assert(graph_node_offsets_.size() == no_graphs());
        assert(graph_no < graph_node_offsets_.size());
        assert(node_no < no_nodes(graph_no));
        return graph_node_offsets_[graph_no] + node_no; 
    }

    std::array<std::size_t,2> multigraph_matching_correlation_clustering_transform::multigraph_matching_from_correlation_clustering_node(const std::size_t i) const
    {
        assert(graph_node_offsets_.size() == no_graphs());
        const std::size_t p = std::lower_bound(graph_node_offsets_.begin(), graph_node_offsets_.end(), i) - graph_node_offsets_.begin();
        assert(p < no_graphs());
        const std::size_t p_node = i - graph_node_offsets_[p];
        assert(i < no_nodes(p));
        assert(i == multigraph_matching_to_correlation_clustering_node(p, p_node));
        return {p, p_node}; 
    }

    std::size_t multigraph_matching_correlation_clustering_transform::compute_no_graphs()
    {
        std::size_t no_graphs = 1 + [](const auto& gm) {
            return std::max(gm.left_graph_no, gm.right_graph_no); 
        }(*std::max_element(mgm->begin(), mgm->end(), [](const auto& gm1, const auto& gm2) { return std::max(gm1.left_graph_no, gm1.right_graph_no) < std::max(gm2.left_graph_no, gm2.right_graph_no); }));

        return no_graphs;
    }

    void multigraph_matching_correlation_clustering_transform::compute_no_nodes()
    {
        no_nodes_.clear();
        no_nodes_.resize(compute_no_graphs(), 0);

        for(auto& gm : *mgm) {
            for(const auto& a : gm.gm_input.assignments) {
                no_nodes_[gm.left_graph_no] = std::max(a.left_node+1, no_nodes_[gm.left_graph_no]);
                no_nodes_[gm.right_graph_no] = std::max(a.right_node+1, no_nodes_[gm.right_graph_no]); 
            }
        }
    }

    void multigraph_matching_correlation_clustering_transform::compute_node_offsets()
    {
        graph_node_offsets_.clear();
        graph_node_offsets_.reserve(no_graphs());
        std::partial_sum(no_nodes_.begin(), no_nodes_.end(), std::back_inserter(graph_node_offsets_));
        std::rotate(graph_node_offsets_.begin(), graph_node_offsets_.begin() + graph_node_offsets_.size()-1, graph_node_offsets_.end());
        assert(graph_node_offsets_[0] == *std::max_element(graph_node_offsets_.begin(), graph_node_offsets_.end()));
        graph_node_offsets_[0] = 0;
        assert(std::is_sorted(graph_node_offsets_.begin(), graph_node_offsets_.end())); 
    }

    void multigraph_matching_correlation_clustering_transform::transform_multigraph_matching_to_correlation_clustering()
    {
        std::vector<std::vector<char>> edge_present_;
        std::vector<double> not_assigned_offset_left;
        std::vector<double> not_assigned_offset_right;

        auto edge_present = [&](const std::size_t matching_no, const std::size_t i, const std::size_t j) -> bool {
            assert(matching_no < mgm->size());
            const auto& gm = (*mgm)[matching_no];
            return edge_present_[matching_no][i * no_nodes(gm.right_graph_no) + j] == 1;
        };

        auto add_edge = [&](const std::size_t matching_no, const std::size_t i, const std::size_t j) -> void {
            assert(matching_no < mgm->size());
            const auto& gm = (*mgm)[matching_no];
            edge_present_[matching_no][i * no_nodes(gm.right_graph_no) + j] = 1;
        };

        auto add_matching = [&](const std::size_t matching_no) -> void {
            assert(matching_no < mgm->size());
            const auto& gm = (*mgm)[matching_no];
            assert(matching_no == edge_present_.size());
            edge_present_.push_back(std::vector<char>(no_nodes(gm.right_graph_no) * no_nodes(gm.left_graph_no),0));;
        };

        // first contruct matching edges
        for(std::size_t matching_no=0; matching_no<mgm->size(); ++matching_no) {
            const auto& gm = (*mgm)[matching_no];
            if(gm.gm_input.quadratic_terms.size() > 0)
                throw std::runtime_error("multigraph matching to multicut transform: can only transform linear multigraph matching problems");

            not_assigned_offset_left.clear();
            not_assigned_offset_left.resize(no_nodes(gm.left_graph_no),0);
            not_assigned_offset_right.clear();
            not_assigned_offset_right.resize(no_nodes(gm.right_graph_no),0);

            assert(gm.left_graph_no < gm.right_graph_no);
            add_matching(matching_no);
            for(const auto& a : gm.gm_input.assignments) {
                if(a.left_node == graph_matching_input::no_assignment) {
                    not_assigned_offset_right[a.right_node] += a.cost;
                } else if(a.right_node == graph_matching_input::no_assignment) {
                    not_assigned_offset_left[a.left_node] += a.cost;
                }
            }
            for(const auto& a : gm.gm_input.assignments) {
                if(a.left_node != graph_matching_input::no_assignment && a.right_node != graph_matching_input::no_assignment) {
                    const std::size_t i = multigraph_matching_to_correlation_clustering_node(gm.left_graph_no, a.left_node);
                    const std::size_t j = multigraph_matching_to_correlation_clustering_node(gm.right_graph_no, a.right_node);
                    cc->add_edge(i, j, a.cost - not_assigned_offset_left[a.left_node] - not_assigned_offset_right[a.right_node]); 
                    add_edge(matching_no, a.left_node, a.right_node);
                }
            }
        }

        const double max_matching_cost = std::abs( std::max_element(cc->begin(), cc->end(), [](const auto e1, const auto& e2) { return std::abs(e1.cost) < std::abs(e2.cost); })->cost );
        const double non_matching_penalty = std::max(1e10, 1000*max_matching_cost);
        for(const auto& a : cc->edges()) { assert(max_matching_cost >= std::abs(a.cost)); }

        // construct negative edges that prevent nodes in one point set to be matched to each other
        for(std::size_t p=0; p<no_graphs(); ++p) {
            for(std::size_t j=0; j<no_nodes(p); ++j) {
                for(std::size_t k=j+1; k<no_nodes(p); ++k) {
                    // TODO: can be made more efficient
                    const std::size_t n1 = multigraph_matching_to_correlation_clustering_node(p,j);
                    const std::size_t n2 = multigraph_matching_to_correlation_clustering_node(p,k);
                    cc->add_edge(n1, n2, non_matching_penalty);
                }
            }
        }

        // if graph matching problem is not dense, add edges that disallow matchings for edges that are not present
        for(std::size_t matching_no=0; matching_no<mgm->size(); ++matching_no) {
            const auto& gm = (*mgm)[matching_no];
            for(std::size_t i=0; i<no_nodes(gm.left_graph_no); ++i) {
                for(std::size_t j=0; j<no_nodes(gm.right_graph_no); ++j) {
                    if(!edge_present(matching_no, i, j)) {
                        const std::size_t n1 = multigraph_matching_to_correlation_clustering_node(gm.left_graph_no,i);
                        const std::size_t n2 = multigraph_matching_to_correlation_clustering_node(gm.right_graph_no,j);
                        cc->add_edge(n1, n2, non_matching_penalty);
                    }
                }
            }
        }
    }

    multigraph_matching_input::labeling multigraph_matching_correlation_clustering_transform::transform(const correlation_clustering_edge_labeling& l_input) const
    {
        multigraph_matching_input::labeling l_output;
        l_output.reserve(mgm->size());

        // TODO: transform cc edge labeling to node labeling and then transform that to multigraph matching labeling
        std::size_t e = 0;
        for(auto& gm : *mgm) {
            linear_assignment_problem_input::labeling l(no_nodes(gm.left_graph_no), std::numeric_limits<std::size_t>::max());
            l_output.push_back( {gm.left_graph_no, gm.right_graph_no, l} );
            for(const auto& a : gm.gm_input.assignments) {
                if(a.left_node != graph_matching_input::no_assignment && a.right_node != graph_matching_input::no_assignment){ 
                    const std::size_t i = multigraph_matching_to_correlation_clustering_node(gm.left_graph_no, a.left_node);
                    const std::size_t j = multigraph_matching_to_correlation_clustering_node(gm.right_graph_no, a.right_node);
                    if(l_input[e] == 1)
                        l_output.back().labeling[a.left_node] = a.right_node;
                    ++e;
                }
            }
        }

        // following edges should disallow infeasible pairwise matchings: they must all be zero.
        for(; e<l_input.size(); ++e)
            if(l_input[e] == 1)
                throw std::runtime_error("correlation clustering not corresponding to multigraph matching.");

        assert(e == l_input.size());
        assert(l_output.check_primal_consistency());

        return l_output;
    }

    correlation_clustering_edge_labeling multigraph_matching_correlation_clustering_transform::transform(const multigraph_matching_input::labeling& input) const
    {
        correlation_clustering_edge_labeling output;

        // first compute a clustering representation
        union_find uf(total_no_nodes());
        for(const auto& gm : input) {
            for(std::size_t i=0; i<gm.labeling.size(); ++i) {
                if(gm.labeling[i] != std::numeric_limits<std::size_t>::max()) {
                    const std::size_t cc_node_1 = multigraph_matching_to_correlation_clustering_node(gm.left_graph_no, i);
                    const std::size_t cc_node_2 = multigraph_matching_to_correlation_clustering_node(gm.right_graph_no, gm.labeling[i]);
                    uf.merge(cc_node_1,cc_node_2);
                }
            }
        }

        for(const auto& e : cc->edges()) {
            const std::size_t i = e[0];
            const std::size_t j = e[1];
            if(uf.connected(i,j))
                output.push_back(1);
            else
                output.push_back(0); 
        }

        return output;
    }

} // namespace LPMP

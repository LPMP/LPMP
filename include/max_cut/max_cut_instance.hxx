#pragma once

#include "cut_base/cut_base_instance.hxx"
#include "max_cut_factors_messages.h"
#include <cassert>
#include "union_find.hxx"

namespace LPMP {
   
// Instead of ordinary max-cut, which is maximized, we take the minimization format to be consistent with max_cut and the LPMP framework

struct max_cut_instance : public cut_base_instance {
    template<typename STREAM>
    void write_problem(STREAM& s) const;
};

struct max_cut_node_labeling : public cut_base_node_labeling {
    using cut_base_node_labeling::cut_base_node_labeling;
};

struct max_cut_edge_labeling : public cut_base_edge_labeling {
    using cut_base_edge_labeling::cut_base_edge_labeling;

    bool check_primal_consistency(const max_cut_instance& instance) const;
    max_cut_node_labeling transform_to_node_labeling(const max_cut_instance& instance) const;
};

struct triplet_max_cut_instance : public triplet_cut_base_instance<max_cut_triplet_factor> {};

struct quadruplet_max_cut_instance : public quadruplet_cut_base_instance<max_cut_quadruplet_factor, triplet_max_cut_instance> {};

inline bool max_cut_edge_labeling::check_primal_consistency(const max_cut_instance& instance) const
{
    assert(this->size() == instance.no_edges());
    // check if edge labeling is valid cut
    union_find uf(instance.no_nodes());
    for(std::size_t e=0; e<instance.edges().size(); ++e)
        if((*this)[e] == 0)
            uf.merge(instance.edges()[e][0], instance.edges()[e][1]);

    for(std::size_t e=0; e<instance.edges().size(); ++e)
        if((*this)[e] == 1)
            if(uf.connected(instance.edges()[e][0], instance.edges()[e][1]))
                return false;

    // check if graph is two-colorable
    graph<char> g(instance.edges().begin(), instance.edges().end(), 
            [&](const auto& e) { 
            const std::size_t pos = std::distance(&instance.edges()[0], &e);
            assert(pos < instance.no_edges());
            return (*this)[pos];
            });

    max_cut_node_labeling output(instance.no_nodes(), 0);
    static constexpr char color_1 = 1;
    static constexpr char color_2 = 2;
    std::vector<char> visited(instance.no_nodes(), 0);
    std::stack<std::size_t> q;
    for(std::size_t n=0; n<visited.size(); ++n) {
        visited[n] = color_1;
        if(!visited[n]) {
            q.push(n);
            while(!q.empty()) {
                const std::size_t i = q.top();
                q.pop();
                for(auto edge_it=g.begin(i); edge_it!=g.end(i); ++edge_it) {
                    const std::size_t j = edge_it->head();
                    const bool cut = edge_it->edge() == 1;
                    if(!visited[j]) {
                        const char color = visited[i];
                        visited[j] = 1;
                        q.push(j);
                        output[j] =  cut ? 1 - output[i] : output[i];
                    } else {
                        if(cut && output[i] == output[j])
                            return false;
                        if(!cut && output[i] == output[j])
                            return false;
                    }
                }
            }
        }
    }

    return true;
}

template<typename STREAM>
void max_cut_instance::write_problem(STREAM& s) const
{
    s << "MAX-CUT\n";
    for(const auto& e : this->edges())
        s << e[0] << " " << e[1] << " " << -e.cost << "\n";
}

inline max_cut_node_labeling max_cut_edge_labeling::transform_to_node_labeling(const max_cut_instance& instance) const
{
    graph<char> g(instance.edges().begin(), instance.edges().end(), 
            [&](const auto& e) { 
            const std::size_t pos = std::distance(&instance.edges()[0], &e);
            assert(pos < instance.no_edges());
            assert( (*this)[pos] == 0 || (*this)[pos] == 1 );
            return (*this)[pos];
            });

    max_cut_node_labeling output(instance.no_nodes(), 0);
    std::vector<char> visited(instance.no_nodes(), 0);
    std::stack<std::size_t> q;
    for(std::size_t n=0; n<visited.size(); ++n) {
        if(!visited[n]) {
            q.push(n);
            while(!q.empty()) {
                const std::size_t i = q.top();
                q.pop();
                for(auto edge_it=g.begin(i); edge_it!=g.end(i); ++edge_it) {
                    const std::size_t j = edge_it->head();
                    if(!visited[j]) {
                        visited[j] = 1;
                        q.push(j);
                        output[j] = edge_it->edge() == 1 ? 1 - output[i] : output[i];
                    }
                }
            }
        }
    }
    return output;
}

} // namespace LPMP

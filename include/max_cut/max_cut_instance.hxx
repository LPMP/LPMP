#ifndef LPMP_MAX_CUT_INSTANCE_HXX
#define LPMP_MAX_CUT_INSTANCE_HXX

#include <vector>
#include <array>
#include <limits>
#include <algorithm>
#include "hash_functions.hxx"
#include "union_find.hxx"

namespace LPMP {

// we take the minimization format
struct max_cut_instance {
    struct weighted_edge : public std::array<std::size_t,2> { 
        weighted_edge(const std::size_t i, const std::size_t j, const double _cost) : std::array<std::size_t,2>({i,j}), cost(_cost) {}
        double cost; 
    };

    double constant = 0.0;
    std::vector<weighted_edge> edges; 

    using labeling = std::vector<unsigned char>;

    double evaluate(const labeling& l) const
    {
        assert(l.size() == edges.size());
        double cost = constant;
        std::size_t no_nodes = 0;
        for(std::size_t e=0; e<edges.size(); ++e) {
            assert(l[e] == 0 || l[e] == 1);
            if(l[e] == 1) {
                cost += edges[e].cost;
            }
            no_nodes = std::max({no_nodes, edges[e][0]+1, edges[e][1]+1});
        }

        // check if edge labeling is valid cut
        union_find uf(no_nodes);
        for(std::size_t e=0; e<edges.size(); ++e) {
            if(l[e] == 0) {
                uf.merge(edges[e][0], edges[e][1]);
            }
        }

        for(std::size_t e=0; e<edges.size(); ++e) {
            if(l[e] == 1) {
                if(uf.connected(edges[e][0], edges[e][1])) {
                    return std::numeric_limits<double>::infinity();
                }
            }
        }

        return cost;
    }

// write out as maximization formulation
    template<typename STREAM>
    void write(STREAM& s) const
    {
        for(const auto& e : edges) {
            s << e[0] << " " << e[1] << " " << -e.cost << "\n";
        }
    }
};

} // namespace LPMP

#endif // LPMP_MAX_CUT_INSTANCE_HXX

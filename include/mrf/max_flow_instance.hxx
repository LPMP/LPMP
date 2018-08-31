#ifndef LPMP_MAX_FLOW_INSTANCE_HXX
#define LPMP_MAX_FLOW_INSTANCE_HXX

#include <vector>
#include <array>
#include <limits>
#include <algorithm>
#include <cassert>

namespace LPMP {

struct max_flow_instance {
    struct capacitated_arc : public std::array<std::size_t,2> { 
        capacitated_arc() {}
        capacitated_arc(const std::size_t i, const std::size_t j, const double c)
        : std::array<std::size_t,2>({i,j}), capacity(c) {}

        double capacity; 
    };

    std::size_t no_nodes;
    std::size_t source, terminal;
    std::vector<capacitated_arc> arcs;

    void add_arc(const std::size_t i, const std::size_t j, const double capacity)
    {
        assert(capacity >= 0.0);
        no_nodes = std::max({no_nodes,i+1,j+1});
        arcs.push_back(capacitated_arc(i,j,capacity));
    }

    using flow = std::vector<long int>;

    double evaluate(const flow& f) const
    {
        assert(f.size() == arcs.size());
        double total_flow = 0.0;
        std::size_t no_nodes = 0;
        for(const auto& a : arcs) {
            if(a[0] == source) {
                total_flow += a.capacity;
            }
            no_nodes = std::max({no_nodes, a[0], a[1]});
        }

        // check if flow is feasible
        std::vector<long int> flow_conservation(no_nodes, 0);
        for(const auto& a : arcs) {
            flow_conservation[a[0]] -= a.capacity;
            flow_conservation[a[1]] += a.capacity;
        }
        for(std::size_t i=0; i<flow_conservation.size(); ++i) {
            if(flow_conservation[i] != 0 && i != source && i != terminal) {
                return std::numeric_limits<double>::infinity(); 
            }
        }
        return total_flow;
    }
};

} // namespace LPMP

#endif // LPMP_MAX_FLOW_INSTANCE_HXX

#ifndef LPMP_MAXFLOW_FACTOR_HXX
#define LPMP_MAXFLOW_FACTOR_HXX

#include "max_flow_instance.hxx"
#include "binary_MRF_instance.hxx"
#include "maxflow/maxflow.h"

namespace LPMP {

class maxflow_factor {
public:

    maxflow_factor(const max_flow_instance& input)
    : m(input.no_nodes-2, input.arcs.size())
    {
        m.add_node(input.no_nodes-2);
        for(const auto& a : input.arcs) {
            assert(a[1] != input.source);
            assert(a[0] != input.terminal);
            if(a[0] == input.source) {
                m.add_tweights(a[1]-2, a.capacity, 0.0); 
            } else if(a[1] == input.terminal) {
                m.add_tweights(a[0]-2, 0.0, a.capacity); 
            } else {
                m.add_edge(a[0]-2,a[1]-2, a.capacity, 0.0);
            } 
        }
    }

    maxflow_factor(const binary_Potts_instance& input)
    : m(input.unaries.size(), input.pairwise_potentials.size())
    {
        assert(input.constant == 0.0);

        m.add_node(input.unaries.size());
        for(std::size_t i=0; i<input.unaries.size(); ++i) {
            m.add_tweights(i, input.unaries[i][1], input.unaries[i][0]); // note: this is correct so, since first comes arc from source to given node, then arc from given node to terminal
        }

        for(const auto& p : input.pairwise_potentials) {
            if(p.cost < 0.0) throw std::runtime_error("graph cut assumes submodular pairwise potentials.");
            m.add_edge(p[0], p[1], p.cost, p.cost);
        }
    }

    double LowerBound() const
    {
        return m.maxflow();
    }

private:
    mutable maxflow::Graph<double,double,double> m;
};

} // namespace LPMP

#endif // LPMP_MAXFLOW_FACTOR_HXX


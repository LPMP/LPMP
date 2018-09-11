#ifndef LPMP_QPBO_FACTOR_HXX
#define LPMP_QPBO_FACTOR_HXX

#include "mrf_input.h"
#include "max_flow_instance.hxx"
#include "binary_MRF_instance.hxx"
#include "qpbo/qpbo.hpp"

namespace LPMP {

class qpbo_factor {
public:

    qpbo_factor(const binary_MRF_instance& input)
    : q(input.unaries.size(), input.pairwise_potentials.size())
    {
        assert(input.constant == 0.0);

        q.AddNode(input.unaries.size());
        for(std::size_t i=0; i<input.unaries.size(); ++i) {
            q.AddUnaryTerm(i, input.unaries[i][0], input.unaries[i][1]);
        }

        for(const auto p : input.pairwise_potentials) {
            q.AddPairwiseTerm(p.i, p.j, p.cost[0][0], p.cost[0][1], p.cost[1][0], p.cost[1][1]);
        }
    }

    qpbo_factor(const binary_Potts_instance& input)
    : q(input.unaries.size(), input.pairwise_potentials.size())
    {
        assert(input.constant == 0.0);

        q.AddNode(input.unaries.size());
        for(std::size_t i=0; i<input.unaries.size(); ++i) {
            q.AddUnaryTerm(i, input.unaries[i][0], input.unaries[i][1]);
        }

        for(const auto p : input.pairwise_potentials) {
            q.AddPairwiseTerm(p[0], p[1], 0.0, p.cost, p.cost, 0.0);
        }
    }

    qpbo_factor(const mrf_input& input)
    : q(input.no_variables(), input.no_pairwise_factors())
    {
       q.AddNode(input.no_variables());
       for(std::size_t i=0; i<input.no_variables(); ++i) {
          if(input.cardinality(i) != 2) throw std::runtime_error("qpbo only accepts binary models.");
          q.AddUnaryTerm(i, input.get_unary(i)[0], input.get_unary(i)[1]);
       }

       for(std::size_t p=0; p<input.no_pairwise_factors(); ++p) {
          const auto vars = input.get_pairwise_variables(p);
          const auto pot = input.get_pairwise_potential(p);
          q.AddPairwiseTerm( vars[0], vars[1], pot(0,0), pot(0,1), pot(1,0), pot(1,1));
       } 
    }

    double LowerBound() const
    {
        q.Solve();
        return q.ComputeTwiceLowerBound()/2.0; 
    }

private:
    mutable qpbo::QPBO<double> q;

};

} // namespace LPMP

#endif // LPMP_QPBO_FACTOR_HXX

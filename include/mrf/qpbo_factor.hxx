#pragma once

#include "mrf_input.h"
#include "max_flow_instance.hxx"
#include "binary_MRF_instance.hxx"
#include "qpbo/qpbo.hpp"
#include <iostream> // TODO: remove

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

    std::size_t no_variables() const { return q.GetNodeNum(); }
    std::size_t no_edges() const { return q.GetMaxEdgeNum(); }

    std::size_t no_persistent_labels() const
    {
       q.Solve();
       q.ComputeWeakPersistencies();
       std::size_t no_persistent_labels = 0;
       for(std::size_t i=0; i<q.GetNodeNum(); ++i)
          if(q.GetLabel(i) >= 0)
             ++no_persistent_labels;
       return no_persistent_labels;
    }

    binary_MRF_instance get_reduced_problem() const
    {
        std::cout << "LB = " << LowerBound() << ", energy = " << q.ComputeTwiceEnergy()/2.0 << "\n";
       binary_MRF_instance output;
       q.Solve();
       output.constant = q.ComputeTwiceLowerBound()/2.0;
       //std::cout << "constant lb model = " << output.lower_bound() << "\n";
       q.ComputeWeakPersistencies();

       std::size_t no_reduced_nodes = 0;
       std::vector<std::size_t> original_to_reduced_node;
       original_to_reduced_node.reserve(no_variables());
       constexpr static std::size_t node_fixed_val = std::numeric_limits<std::size_t>::max();
       auto node_fixed = [&](const std::size_t i) {
          assert(i < original_to_reduced_node.size());
          return original_to_reduced_node[i] == node_fixed_val; 
       };
       auto get_reduced_node = [&](const std::size_t i) { 
          assert(!node_fixed(i));
          return original_to_reduced_node[i]; 
       };

       for(std::size_t i=0; i<q.GetNodeNum(); ++i) {
           double E0, E1;
           q.GetTwiceUnaryTerm(i, E0, E1);
           // substraction by below value keeps lower bound intact
           E0 = E0 - std::min(E0,E1);
           E1 = E1 - std::min(E0,E1);
           if(q.GetLabel(i) < 0) { // node persistently labelled
               output.unaries.push_back({0.5*E0, 0.5*E1});
               original_to_reduced_node.push_back(no_reduced_nodes++);
           } else {
               assert(q.GetLabel(i) == 0 || q.GetLabel(i) == 1);
               original_to_reduced_node.push_back(node_fixed_val);
               output.constant += 0.5 * (q.GetLabel(i) == 0 ? E0 : E1);
           }
       }

       for(std::size_t e=0; e<no_edges(); ++e) {
           int i, j;
           std::array<std::array<double,2>,2> E;
           q.GetTwicePairwiseTerm(e, i, j, E[0][0], E[0][1], E[1][0], E[1][1]);
           binary_MRF_instance::binary_pairwise_potential p(i,j, E);
           p.cost[0][0] *= 0.5;
           p.cost[0][1] *= 0.5;
           p.cost[1][0] *= 0.5;
           p.cost[1][1] *= 0.5;
           if(!node_fixed(p.i) && !node_fixed(p.j)) {
               p.i = get_reduced_node(p.i);
               p.j = get_reduced_node(p.j);
               output.pairwise_potentials.push_back(p);
           } else if(!node_fixed(p.i) && node_fixed(p.j)) {
               const std::size_t reduced_i = get_reduced_node(p.i);
               output.unaries[reduced_i][0] += p.cost[0][q.GetLabel(p.j)];
               output.unaries[reduced_i][1] += p.cost[1][q.GetLabel(p.j)];
           } else if(node_fixed(p.i) && !node_fixed(p.j)) {
               const std::size_t reduced_j = get_reduced_node(p.j);
               output.unaries[reduced_j][0] += p.cost[q.GetLabel(p.i)][0];
               output.unaries[reduced_j][1] += p.cost[q.GetLabel(p.i)][1];
           } else {
               output.constant += p.cost[q.GetLabel(p.i)][q.GetLabel(p.j)]; 
           }
       }

       output.merge_parallel_edges();
       return output;
    }

private:
    mutable qpbo::QPBO<double> q;

};

} // namespace LPMP

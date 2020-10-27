#pragma once

#include <cassert>
#include "mrf/binary_MRF_instance.hxx"
#include "max_cut/max_cut_instance.hxx"

namespace LPMP {

max_cut_instance transform_binary_Potts_to_max_cut(const binary_Potts_instance& input)
{
    max_cut_instance output;
    output.edges().reserve(input.pairwise_potentials.size() + input.unaries.size());

    for(const auto& p : input.pairwise_potentials) {
        output.edges().push_back({p[0]+1, p[1]+1, p.cost});
    }

    for(std::size_t i=0; i<input.unaries.size(); ++i) {
        const auto& u = input.unaries[i];
        output.edges().push_back({0, i+1, u[1] - u[0]});
        output.add_to_constant(u[0]);
    }

    return output; 
}

max_cut_edge_labeling transform_binary_Potts_labeling_to_max_cut(const binary_Potts_instance& instance, const binary_Potts_instance::labeling& input)
{
    max_cut_edge_labeling l;
    l.reserve(instance.unaries.size() + instance.pairwise_potentials.size());

    assert(input.size() == instance.unaries.size());
    for(const auto l : input) { assert(l == 0 || l == 1); }

    for(const auto& p : instance.pairwise_potentials) {
        const std::size_t i = p[0];
        const std::size_t j = p[1];
        const auto label_i = input[i];
        const auto label_j = input[j];
        const unsigned char label = (label_i == label_j) ? 0 : 1;
        l.push_back(label);
    }

    for(std::size_t i=0; i<instance.unaries.size(); ++i) {
        l.push_back(input[i]);
    }

    return l;
}

} // namespace LPMP

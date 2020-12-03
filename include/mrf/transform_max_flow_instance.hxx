#ifndef LPMP_TRANSFORM_MAX_FLOW_HXX
#define LPMP_TRANSFORM_MAX_FLOW_HXX

#include "max_flow_instance.hxx"
#include "binary_MRF_instance.hxx"

namespace LPMP {

// every node in MRF corresponds to two consecutive nodes (first even, second odd) in the max flow graph.
binary_MRF_instance transform_QPBO_max_flow_to_binary_MRF(const max_flow_instance& instance)
{
    assert(instance.source <= 1 && instance.terminal <= 1 && instance.source != instance.terminal);
    assert(instance.no_nodes >= 2 && instance.no_nodes%2 == 0);

    binary_MRF_instance binary_MRF;
    binary_MRF.unaries.resize((instance.no_nodes-1)/2, {0.0,0.0});

    // returns (node, label)
    auto get_mrf_node_label = [&instance](const std::size_t max_flow_node) -> std::tuple<std::size_t, std::size_t> { 
       if(max_flow_node <= 1) throw std::runtime_error("max flow node does not correspond to mrf node, but is source or sink.");
       assert(max_flow_node > 1 && max_flow_node < instance.no_nodes);
       const std::size_t no_mrf_nodes = (instance.no_nodes-2)/2;
       const std::size_t mrf_node = (max_flow_node - 2) % no_mrf_nodes;
       const std::size_t label = max_flow_node - 2 < no_mrf_nodes ? 0 : 1;
       return {mrf_node, label};
    };

    // first come the unaries
    std::size_t i=0;
    for(; i<instance.arcs.size(); ++i) {
        const auto& arc = instance.arcs[i];
        const auto cap = arc.capacity;

        if(arc[0] == instance.source) {
            if(arc[1] == instance.terminal) throw std::runtime_error("arc from source to terminal not allowed.");
            assert(arc[1] >= 2);
            const auto [mrf_node, label] = get_mrf_node_label(arc[1]);
            assert(mrf_node < binary_MRF.unaries.size());
            binary_MRF.unaries[mrf_node][1-label] += cap;
        } else if(arc[1] == instance.terminal) {
            if(arc[0] == instance.source) throw std::runtime_error("arc from source to terminal not allowed.");
            assert(arc[0] >= 2);
            const auto [mrf_node, label] = get_mrf_node_label(arc[0]);
            assert(mrf_node < binary_MRF.unaries.size());
            binary_MRF.unaries[mrf_node][label] += cap; 
        } else {
            break;
        }
    }   

    // now pairwise potentials begin
    // pairwise potentials come in pairs of two sister edges (but not necessarily sequentially).
    if((instance.arcs.size() - i)%4 != 0) throw std::runtime_error("pairwise potential arcs must come in groups of four.");

    std::vector<binary_MRF_instance::binary_pairwise_potential> pairwise_potentials;
    pairwise_potentials.reserve((instance.arcs.size()-i)/2);

    for(; i<instance.arcs.size(); ++i) {
        const auto [mrf_node_1, label_1] = get_mrf_node_label(instance.arcs[i][0]);
        const auto [mrf_node_2, label_2] = get_mrf_node_label(instance.arcs[i][1]);
        if(mrf_node_1 == mrf_node_2) throw std::runtime_error("max flow arc between identical mrf node.");

        binary_MRF_instance::binary_pairwise_potential p(mrf_node_1, mrf_node_2, {{{0.0,0.0},{0.0,0.0}}} );
        auto& pot = p.cost;

        const auto cap = instance.arcs[i].capacity;
        assert(cap >= 0.0);

        if(label_1 == 0 && label_2 == 0) {
           assert(pot[0][1] == 0.0);
           pot[0][1] += cap;
        } else if(label_1 == 0 && label_2 == 1) {
           assert(pot[0][0] == 0.0);
           pot[0][0] += cap;
        } else if(label_1 == 1 && label_2 == 0) {
           assert(pot[1][1] == 0.0);
           pot[1][1] += cap;
        } else if(label_1 == 1 && label_2 == 1) { 
           assert(pot[1][0] == 0.0);
           pot[1][0] += cap;
        } else {
           assert(false); 
        }

        if(mrf_node_2 < mrf_node_1) {
           p.transpose();
        }
        pairwise_potentials.push_back(p);
    }

    // merge pairwise potentials on same variables
    std::sort(pairwise_potentials.begin(), pairwise_potentials.end());
    binary_MRF.pairwise_potentials.reserve(pairwise_potentials.size()/2);

    assert( (pairwise_potentials.size()%4) == 0);

    for(std::size_t i=0; i<pairwise_potentials.size(); i+=4) {
       if( !pairwise_potentials[i].has_same_support(pairwise_potentials[i+1])
             || !pairwise_potentials[i+1].has_same_support(pairwise_potentials[i+2])
             || !pairwise_potentials[i+2].has_same_support(pairwise_potentials[i+3]) ) 
       {
          throw std::runtime_error("two pairs of sister arcs that either cross (nonsubmodular) or are parallel (submodular) must be present.");
       }
       auto pot = pairwise_potentials[i];
       pot.add_cost(pairwise_potentials[i+1]);
       pot.add_cost(pairwise_potentials[i+2]);
       pot.add_cost(pairwise_potentials[i+3]);
       binary_MRF.pairwise_potentials.push_back(pot);
    }

    // check that there are no two potentials on the same variables.
    std::sort(binary_MRF.pairwise_potentials.begin(), binary_MRF.pairwise_potentials.end());
    if(binary_MRF.pairwise_potentials.size() > 0) {
       for(std::size_t i=0; i<binary_MRF.pairwise_potentials.size()-1; ++i) {
          const auto& p1 = binary_MRF.pairwise_potentials[i];
          const auto& p2 = binary_MRF.pairwise_potentials[i+1];
          if(p1.has_same_support(p2)) throw std::runtime_error("there must not be more than one pairwise potentials for every combination of variables."); 
       }
    }

    return binary_MRF;
}

binary_Potts_instance transform_graph_cut_max_flow_to_binary_Potts(const max_flow_instance& input)
{
    assert(input.source <= 1);
    assert(input.terminal <= 1);
    assert(input.no_nodes >= 2);

    binary_Potts_instance binary_MRF;
    binary_MRF.unaries.resize(input.no_nodes-2, {0.0, 0.0});

    // first come the unaries
    std::size_t i=0;
    for(; i<input.arcs.size(); ++i) {
        const auto& arc = input.arcs[i];
        if(arc[0] == input.source) {
            assert(arc[1] >= 2);
            const auto node = arc[1]-2;
            assert(node < binary_MRF.unaries.size());
            binary_MRF.unaries[node][0] += arc.capacity;
        } else if(arc[1] == input.terminal) {
            assert(arc[0] >= 2);
            const auto node = arc[0]-2;
            assert(node < binary_MRF.unaries.size());
            binary_MRF.unaries[node][1] += arc.capacity;
        } else {
            break;
        }
    }

    // now pairwise potentials begin in group of two arcs
    assert((input.arcs.size() - i)%2 == 0);

    for(; i<input.arcs.size(); i+=2) {
        const std::size_t node_1 = std::min({ input.arcs[i][0], input.arcs[i][1], input.arcs[i+1][0], input.arcs[i+1][1] });
        const std::size_t node_2 = std::max({ input.arcs[i][0], input.arcs[i][1], input.arcs[i+1][0], input.arcs[i+1][1] });
        assert(node_1 > 1);
        assert(node_1 != node_2);

        binary_MRF_instance::binary_pairwise_potential p;

        for(std::size_t j=0; j<2; ++j) {
            const std::size_t u=input.arcs[i+j][0];
            const std::size_t v=input.arcs[i+j][1];
            assert(u == node_1 || u == node_2);
            assert(v == node_1 || v == node_2);
            const auto cost = input.arcs[i+j].capacity;
            if(u == node_1) {
                p.cost[1][0] += cost; 
            } else  {
                assert(u == node_2);
                p.cost[0][1] += cost; 
            }
        } 
        auto [c, msg_1, msg_2] = p.make_potts();
        assert(c >= 0.0); // must be submodular

        binary_MRF.unaries[node_1-2][0] += msg_1[0];
        binary_MRF.unaries[node_1-2][1] += msg_1[1];

        binary_MRF.unaries[node_2-2][0] += msg_2[0];
        binary_MRF.unaries[node_2-2][1] += msg_2[1];

        binary_Potts_instance::weighted_edge e({node_1-2,node_2-2,c}); 
        binary_MRF.pairwise_potentials.push_back(e);
    }

    return binary_MRF;
}

} // namespace LPMP

#endif // LPMP_TRANSFORM_MAX_FLOW_HXX

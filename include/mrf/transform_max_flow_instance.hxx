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

    auto get_mrf_node = [](const std::size_t max_flow_node) { assert(max_flow_node > 1); return max_flow_node/2 - 1; };

    // first come the unaries
    std::size_t i=0;
    for(; i<instance.arcs.size(); ++i) {
        const auto& arc = instance.arcs[i];
        const auto cap = arc.capacity;

        if(arc[0] == instance.source) {
            if(arc[1] == instance.terminal) throw std::runtime_error(" arc from source to terminal not allowed.");
            assert(arc[1] >= 2);
            const std::size_t node = get_mrf_node(arc[1]);
            assert(node < binary_MRF.unaries.size());
            if(arc[1] % 2 == 0) {
                binary_MRF.unaries[node][0] += cap;
            } else {
                binary_MRF.unaries[node][1] += cap;
            }
        } else if(arc[1] == instance.terminal) {
            if(arc[0] == instance.source) throw std::runtime_error(" arc from source to terminal not allowed.");
            assert(arc[0] >= 2);
            const auto node = get_mrf_node(arc[0]);
            assert(node < binary_MRF.unaries.size());
            if(arc[0] % 2 == 0) {
                binary_MRF.unaries[node][1] += cap;
            } else {
                binary_MRF.unaries[node][0] += cap;
            }
        } else {
            break;
        }
    }   

    // now pairwise potentials begin in groups of four arcs
    if((instance.arcs.size() - i)%4 != 0) throw std::runtime_error("pairwise potential arcs must come in groups of four.");

    for(; i<instance.arcs.size(); i+=4) {
        std::array<std::size_t,8> nodes = {
            instance.arcs[i][0], instance.arcs[i][1], 
            instance.arcs[i+1][0], instance.arcs[i+1][1],
            instance.arcs[i+2][0], instance.arcs[i+2][1],
            instance.arcs[i+3][0], instance.arcs[i+3][1] 
            }; 

        const std::size_t node_1 = *std::min_element(nodes.begin(), nodes.end());
        const std::size_t node_2 = *std::max_element(nodes.begin(), nodes.end()) - 1;
        if(node_1 <= 1) throw std::runtime_error("edge from source or to terminal not allowed in pairwise section of max flow input.");
        assert(node_1 < node_2);

        for(std::size_t j=0; j<4; ++j) {
            const std::size_t u=instance.arcs[i+j][0];
            const std::size_t v=instance.arcs[i+j][1];
            assert(u == node_1 || u == node_1+1 || u == node_2 || u == node_2+1);
            assert(v == node_1 || v == node_1+1 || v == node_2 || v == node_2+1);
        } 

        binary_MRF_instance::binary_pairwise_potential p(get_mrf_node(node_1), get_mrf_node(node_2), {{{0.0,0.0},{0.0,0.0}}} );
        auto& pot = p.cost;

        for(std::size_t j=0; j<4; ++j) {
            const std::size_t u = instance.arcs[i+j][0];
            const std::size_t v = instance.arcs[i+j][1];
            const auto cap = instance.arcs[i+j].capacity;
            assert(cap >= 0.0);

            if(u == node_1 && v == node_2) {
                pot[1][0] += cap; 
            } else if(u == node_1 && v == node_2+1) {
                pot[1][1] += cap;
            } else if(u == node_1+1 && v == node_2) {
                pot[0][0] += cap;
            } else if(u == node_1+1 && v == node_2+1) {
                pot[0][1] += cap; 

            } else if(u == node_2 && v == node_1) {
                pot[0][1] += cap;
            } else if(u == node_2 && v == node_1+1) {
                pot[1][1] += cap;
            } else if(u == node_2+1 && v == node_1) {
                pot[0][0] += cap; 
            } else if(u == node_2+1 && v == node_1+1) {
                pot[1][0] += cap; 
            }  else {
                throw std::runtime_error("max flow input must have four consecutive arcs from a pair of involution nodes.");
            }   
        }

        binary_MRF.pairwise_potentials.push_back(p);
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

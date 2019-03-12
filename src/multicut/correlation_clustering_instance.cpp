#include "multicut/correlation_clustering_instance.h"

namespace LPMP {

multicut_instance correlation_clustering_instance::transform_to_multicut() const
{
    multicut_instance output;
    for(auto& e : edges()) {
        output.add_to_constant(e.cost);
        output.add_edge(e[0], e[1], -e.cost[0]);
    } 
    return output;
}
correlation_clustering_edge_labeling correlation_clustering_node_labeling::transform_to_edge_labeling(const correlation_clustering_instance& instance) const 
{
   correlation_clustering_edge_labeling output;
   output.reserve(instance.no_edges());
   for(const auto& e : instance.edges()) {
      const bool cut = (*this)[e[0]] != (*this)[e[1]];
      output.push_back(cut); 
   }

   assert(std::abs(instance.evaluate(*this) - instance.evaluate(output)) <= 1e-8);
   return output;
} 

bool correlation_clustering_edge_labeling::check_primal_consistency(const correlation_clustering_instance& input) const
{
   union_find uf(input.no_nodes());
   for(std::size_t e=0; e<this->size(); ++e)
      if((*this)[e] == 0)
         uf.merge(input.edges()[e][0], input.edges()[e][1]);
   for(std::size_t e=0; e<this->size(); ++e)
      if((*this)[e] == 0 && uf.find(input.edges()[e][0]) == uf.find(input.edges()[e][1]))
         return false;
   return true;
}

correlation_clustering_node_labeling correlation_clustering_edge_labeling::transform_to_node_labeling(const correlation_clustering_instance& instance) const
{
    assert(this->size() == instance.no_edges());

    union_find uf(instance.no_nodes());

    for(std::size_t e=0; e<this->size(); ++e)
        if((*this)[e] == 0)
            uf.merge(instance.edges()[e][0], instance.edges()[e][1]); 

    correlation_clustering_node_labeling output(instance.no_nodes());
    for(std::size_t i=0; i<instance.no_nodes(); ++i)
        output[i] = uf.find(i);

    assert(std::abs(instance.evaluate(*this) - instance.evaluate(output)) <= 1e-8);
    return output;
}

multicut_edge_labeling correlation_clustering_edge_labeling::transform_to_multicut() const
{
    multicut_edge_labeling output;
    output.reserve(this->size());
    for(const auto x : *this)
        output.push_back(1-x);
    return output;
}

} // namespace LPMP

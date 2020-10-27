#include "multicut/multicut_kernighan_lin.h"
//#include "andres/graph/multicut-lifted/BEC.hxx"
//#include "andres/graph/multicut-lifted/BEC_cut.hxx"

namespace LPMP {

   std::pair<andres::graph::Graph<>, std::vector<double>> construct_andres_multicut_instance(const multicut_instance& instance)
   {
      andres::graph::Graph<> graph(instance.no_nodes());
      std::vector<double> edge_values;
      edge_values.reserve(instance.no_edges());
      for(const auto& e : instance.edges()) {
         graph.insertEdge(e[0], e[1]); 
         edge_values.push_back(e.cost); // TODO: is maximization or minimization done?
      }

      return {graph, edge_values};
   }

   multicut_edge_labeling compute_gaec(const multicut_instance& instance)
   {
      auto [graph, edge_values] = construct_andres_multicut_instance(instance);

      multicut_edge_labeling labeling;
      labeling.resize(instance.no_edges(), 1);

      if(graph.numberOfEdges() > 0)
          andres::graph::multicut::greedyAdditiveEdgeContraction(graph, edge_values, labeling);
      //andres::graph::multicut_lifted::balancedEdgeContraction(graph, graph, edge_values, labeling);
      //andres::graph::multicut_lifted::balancedEdgeContraction_cut(graph, graph, edge_values, labeling);

      std::cout << "multicut cost on edges = " << instance.evaluate(labeling) << "\n";

      return labeling;
   }

   multicut_edge_labeling compute_multicut_kernighan_lin(const multicut_instance& instance, multicut_edge_labeling labeling)
   { 
      auto [graph, edge_values] = construct_andres_multicut_instance(instance);

      if(labeling.size() != instance.no_edges()) {
         labeling.clear();
         labeling.resize(instance.no_edges(), 1);
      }

      if(graph.numberOfEdges() > 0)
         andres::graph::multicut::kernighanLin(graph, edge_values, labeling, labeling);

      return labeling; 
   }

   multicut_edge_labeling compute_multicut_gaec_kernighan_lin(const multicut_instance& instance)
   {
      auto [graph, edge_values] = construct_andres_multicut_instance(instance);

      multicut_edge_labeling labeling;
      labeling.resize(graph.numberOfEdges(),1);

      if(graph.numberOfEdges() > 0) {
         andres::graph::multicut::greedyAdditiveEdgeContraction(graph, edge_values, labeling);
         andres::graph::multicut::kernighanLin(graph, edge_values, labeling, labeling);
      }

      return labeling;
   }

   multicut_edge_labeling compute_multicut_greedy_edge_fixation(const multicut_instance& instance)
   {
      auto [graph, edge_values] = construct_andres_multicut_instance(instance);
      multicut_edge_labeling labeling;
      labeling.resize(graph.numberOfEdges(),1);
      andres::graph::multicut::greedyFixation(graph, edge_values, labeling);
      return labeling; 
   }

} // namespace LPMP

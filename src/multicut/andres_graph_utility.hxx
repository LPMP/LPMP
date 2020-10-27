#ifndef GRAPH_ANDRES_UTILITY_HXX
#define GRAPH_ANDRES_UTILITY_HXX

#include <fstream>
#include "andres/graph/graph.hxx"
#include "andres/graph/grid-graph.hxx"
#include "andres/graph/hdf5/graph.hxx"
#include "andres/graph/hdf5/grid-graph.hxx"
#include "andres/functional.hxx"

template<typename GRAPH_TYPE> // may be andres::graph::GridGraph<N> or andres::graph::Graph
void read_graph(const std::string& filename, const std::string& output_file)
{
   GRAPH_TYPE originalGraph;
   andres::graph::Graph<> liftedGraph;
   std::vector<double> edgeValues;

   auto fileHandle = andres::graph::hdf5::openFile(filename);
   andres::graph::hdf5::load(fileHandle, "graph", originalGraph);
   andres::graph::hdf5::load(fileHandle, "graph-lifted", liftedGraph);

   std::vector<size_t> shape;
   andres::graph::hdf5::load(fileHandle, "edge-cut-probabilities", shape, edgeValues);
   andres::graph::hdf5::closeFile(fileHandle);

   // transform to energy cost
   std::transform(edgeValues.begin(), edgeValues.end(), edgeValues.begin(), andres::NegativeLogProbabilityRatio<double,double>());
   assert(edgeValues.size() == liftedGraph.numberOfEdges());

   std::ofstream file_stream(output_file, std::ofstream::out);

   file_stream << "MULTICUT\n";
   for(std::size_t e=0; e<liftedGraph.numberOfEdges(); ++e) {
     auto i = liftedGraph.vertexOfEdge(e,0);
     auto j = liftedGraph.vertexOfEdge(e,1);
     assert(i<j);
     if(originalGraph.findEdge(i,j).first) {
       file_stream << i << " " << j << " " << edgeValues[e] << "\n";
     }
   }

   if(liftedGraph.numberOfEdges() > 0) {
      file_stream << "LIFTED\n"; 
      for(std::size_t e=0; e<liftedGraph.numberOfEdges(); ++e) {
         auto i = liftedGraph.vertexOfEdge(e,0);
         auto j = liftedGraph.vertexOfEdge(e,1);
         assert(i<j);
         if(!originalGraph.findEdge(i,j).first) {
            file_stream << i << " " << j << " " << edgeValues[e] << "\n";
         }
      }
   }
}

#endif

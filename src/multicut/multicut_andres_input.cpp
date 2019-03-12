#include "multicut/multicut_andres_input.h"
#include "andres/graph/graph.hxx"
#include "andres/graph/grid-graph.hxx"
#include "andres/graph/hdf5/graph.hxx"
#include "andres/graph/hdf5/grid-graph.hxx"
#include "andres/functional.hxx"

namespace LPMP {

namespace multicut_andres_input {

   template<typename GRAPH_TYPE> // may be andres::graph::GridGraph<N> or andres::graph::Graph<>
   multicut_instance read_graph(const std::string& filename)
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

      multicut_instance output;

      for(std::size_t e=0; e<liftedGraph.numberOfEdges(); ++e) {
         auto i = liftedGraph.vertexOfEdge(e,0);
         auto j = liftedGraph.vertexOfEdge(e,1);
         assert(i<j);
         if(originalGraph.findEdge(i,j).first) {
            output.add_edge(i, j, edgeValues[e]);
         }
      }

      // TODO: not implemented yet
      assert(liftedGraph.numberOfEdges() == 0);
      /*
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
      */

      return output;
   }


   multicut_instance parse_file(const std::string& filename)
   {
      return read_graph<andres::graph::Graph<>>(filename);
   }

   multicut_instance parse_file_2d_grid_graph(const std::string& filename)
   {
      return read_graph<andres::graph::GridGraph<2>>(filename);
   }

}

}

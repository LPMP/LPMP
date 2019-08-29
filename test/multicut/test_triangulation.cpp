#include "multicut/multicut_instance.h"
#include "multicut/multicut_cycle_packing_parallel.h"
#include "../src/multicut/multicut_cycle_packing_parallel.cpp"
#include "test.h"

using namespace LPMP;


void print_triangle_to_edge(std::vector<triangle_item> triangle_to_edge){
	for (auto t: triangle_to_edge){
         std::cout << "Triangle nodes: " << t.nodes[0] << ", " << t.nodes[1] << ", " << t.nodes[2];
         std::cout << ".  Weights: "  << t.weights[0] << ", " << t.weights[1] << ", " << t.weights[2];
         std::cout << ".  Edge indices: "  << t.edge_indices[0] << ", " << t.edge_indices[1] << ", " << t.edge_indices[2] << std::endl;
      }
}

void print_edge_to_triangle(std::vector<edge_item> edge_to_triangle){
	std::cout << "edge_to_triangle:\n"; 
	for (auto e: edge_to_triangle){
		std::cout << "Edge nodes: " << e.nodes[0] << ", " << e.nodes[1] << ". Cost: " << e.cost;
		for (auto e_t :  e.triangle_indices)
			std::cout << " The third node: " << e_t[0] << ", in triangle: " << e_t[1] << ".  ";
		std::cout << "\n";
     }
}

void print_other_edges(std::vector<edge_t> other_edges){
	std::cout << "Other edges:\n"; 
	for (auto e: other_edges)
		std::cout << e[0] << ", " << e[1] << ": " << e.cost << std::endl; 
}

int main()
{
   {
		multicut_instance test_instance;
		test_instance.add_edge(0,1,1);
		test_instance.add_edge(0,2,1);
		test_instance.add_edge(1,2,-3);
		test_instance.add_edge(1,3,1);
		test_instance.add_edge(2,3,1);
		std::vector<triangle_item> triangle_to_edge;
		std::vector<edge_item> edge_to_triangle;
		std::vector<edge_t> other_edges;

		// check whether triangulation is correct
		multicut_cycle_packing_parallel_impl(test_instance, false, 1, true, triangle_to_edge, edge_to_triangle, other_edges);
		test(triangle_to_edge.size() == 2);
		test(edge_to_triangle.size() == 5);
		print_triangle_to_edge(triangle_to_edge);
		print_edge_to_triangle(edge_to_triangle);
		print_other_edges(other_edges);
   
   }
   {
   		std::cout << "\n==============================\n";
		multicut_instance test_instance;
		test_instance.add_edge(0,1,1);
		test_instance.add_edge(0,2,1);
		test_instance.add_edge(1,2,-1.5);
		test_instance.add_edge(1,3,1);
		test_instance.add_edge(2,3,1);
		std::vector<triangle_item> triangle_to_edge;
		std::vector<edge_item> edge_to_triangle;
		std::vector<edge_t> other_edges;

		// check whether triangulation is correct
		multicut_cycle_packing_parallel_impl(test_instance, false, 1, true, triangle_to_edge, edge_to_triangle, other_edges);
		test(triangle_to_edge.size() == 2);
		test(edge_to_triangle.size() == 5);
		print_triangle_to_edge(triangle_to_edge);
		print_edge_to_triangle(edge_to_triangle);
   		print_other_edges(other_edges);
   }
   {
   		std::cout << "\n==============================\n";
		multicut_instance test_instance;
		test_instance.add_edge(0,1,1);
		test_instance.add_edge(0,2,1);
		test_instance.add_edge(1,2,-1.5);
		test_instance.add_edge(1,3,1);
		test_instance.add_edge(2,3,1);
		test_instance.add_edge(2,4,1);
		std::vector<triangle_item> triangle_to_edge;
		std::vector<edge_item> edge_to_triangle;
		std::vector<edge_t> other_edges;

		// check whether triangulation is correct
		multicut_cycle_packing_parallel_impl(test_instance, false, 1, true, triangle_to_edge, edge_to_triangle, other_edges);
		test(triangle_to_edge.size() == 2);
		test(edge_to_triangle.size() == 5);
		print_triangle_to_edge(triangle_to_edge);
		print_edge_to_triangle(edge_to_triangle);
		print_other_edges(other_edges);
   
   }

}
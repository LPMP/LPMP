#include "test.h"
#include "graph.hxx"
#include <unordered_map>

using namespace LPMP;

int main(int argc, char** argv)
{
	std::vector<std::array<std::size_t,2>> edges = { {0,1}, {1,2}, {2,3}, {0,3}, {0,2} };
	std::sort(edges.begin(), edges.end(), [](const auto& e1, const auto& e2) { return std::lexicographical_compare(e1.begin(), e1.end(), e2.begin(), e2.end()); });

	struct empty {};
	graph<empty> g(edges.begin(), edges.end());

	test(g.no_nodes() == 4);
	test(g.no_edges(0) == 3);
	test(g.no_edges(1) == 2);
	test(g.no_edges(2) == 3);
	test(g.no_edges(3) == 2);

	for(auto edge : edges) {
		test(g.edge_present(edge[0], edge[1]));
	}

	std::vector<std::array<std::size_t,2>> edges_check;
	g.for_each_edge([&](const std::size_t i, const std::size_t j, const empty&) { edges_check.push_back({i,j}); });
	std::sort(edges_check.begin(), edges_check.end(), [](const auto& e1, const auto& e2) { return std::lexicographical_compare(e1.begin(), e1.end(), e2.begin(), e2.end()); });
	test(edges == edges_check); 

	std::vector<std::array<std::size_t,3>> triangles;
	g.for_each_triangle([&](const std::size_t i, const std::size_t j, const std::size_t k, const auto& ij, const auto& ik, const auto& jk) { triangles.push_back({i,j,k}); });
	test(triangles.size() == 2);
	std::sort(triangles.begin(), triangles.end(), [](const auto& t1, const auto& t2) { return std::lexicographical_compare(t1.begin(), t1.end(), t2.begin(), t2.end()); });
	test(triangles[0][0] == 0 && triangles[0][1] == 1 && triangles[0][2] == 2);
	test(triangles[1][0] == 0 && triangles[1][1] == 2 && triangles[1][2] == 3);

	std::vector<std::array<std::size_t,4>> quadrangles;
	g.for_each_quadrangle([&](const std::size_t i, const std::size_t j, const std::size_t k, const std::size_t l) { 
			std::array<std::size_t,4> q({i,j,k,l});	
			std::sort(q.begin(), q.end());
			quadrangles.push_back(q);
			});
	test(quadrangles.size() == 1);
	test(quadrangles[0][0] == 0 && quadrangles[0][1] == 1 && quadrangles[0][2] == 2 && quadrangles[0][3] == 3); 

	decltype(edges) contraction_edges({{0,2}});
	auto [contracted_graph, contraction_mapping] = g.contract(contraction_edges.begin(), contraction_edges.end());
	test(contracted_graph.no_nodes() == 3);
}

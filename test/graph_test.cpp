#include "test.h"
#include "graph.hxx"
#include "dynamic_graph.hxx"
#include "dynamic_graph_thread_safe.hxx"
#include <atomic>

using namespace LPMP;

static std::vector<std::array<std::size_t,2>> edges = { {0,1}, {1,2}, {2,3}, {0,3}, {0,2} };

std::vector<std::vector<size_t>> maximal_cliques = {
    {0,1,2},
    {0,2,3}
};

struct empty {};

template<typename GRAPH_TYPE>
GRAPH_TYPE construct_and_test_graph()
{
	GRAPH_TYPE g(edges.begin(), edges.end());

	test(g.no_nodes() == 4);
	test(g.no_edges(0) == 3); 
	test(g.no_edges(1) == 2);
	test(g.no_edges(2) == 3);
	test(g.no_edges(3) == 2);

    test(g.edge_present(0,1)); test(g.edge_present(1,0));
    test(g.edge_present(1,2)); test(g.edge_present(2,1));
    test(g.edge_present(2,3)); test(g.edge_present(3,2));
    test(g.edge_present(0,3)); test(g.edge_present(3,0));
    test(g.edge_present(0,2)); test(g.edge_present(2,0));

    test(!g.edge_present(1,3)); test(!g.edge_present(3,1));

	for(auto e : edges) {
		test(g.edge_present(e[0], e[1]));
		test(g.edge_present(e[1], e[0]));
	}

    return g;
}

void test_maximal_clique_enumeration()
{
    graph<empty> g(edges.begin(), edges.end());
    std::vector<bool> clique_visited(maximal_cliques.size(), false);
    auto check_clique = [&](auto v_begin, auto v_end) {
        std::vector<size_t> v(v_begin, v_end);
        std::sort(v.begin(), v.end());
        const size_t clique_no = std::distance(maximal_cliques.begin(), std::find(maximal_cliques.begin(), maximal_cliques.end(), v));
        test(v == maximal_cliques[0] || v == maximal_cliques[1]); 
        test(clique_no < maximal_cliques.size());
        test(v == maximal_cliques[clique_no]);
        test(clique_visited[clique_no] == false);
        clique_visited[clique_no] = true; 
    };
    g.for_each_maximal_clique(check_clique);
    test(std::count(clique_visited.begin(), clique_visited.end(), false) == 0);
}

int main(int argc, char** argv)
{
	std::sort(edges.begin(), edges.end(), [](const auto& e1, const auto& e2) { return std::lexicographical_compare(e1.begin(), e1.end(), e2.begin(), e2.end()); });

    construct_and_test_graph<dynamic_graph<empty>>();
    construct_and_test_graph<dynamic_graph_thread_safe<empty>>();
    auto g = construct_and_test_graph<graph<empty>>();

    test_maximal_clique_enumeration();

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
	g.for_each_quadrangle([&](std::array<std::size_t,4> nodes) { 
			std::sort(nodes.begin(), nodes.end());
			quadrangles.push_back(nodes);
			});
	test(quadrangles.size() == 1);
	test(quadrangles[0][0] == 0 && quadrangles[0][1] == 1 && quadrangles[0][2] == 2 && quadrangles[0][3] == 3); 

	decltype(edges) contraction_edges({{0,2}});
	auto [contracted_graph, contraction_mapping] = g.contract(contraction_edges.begin(), contraction_edges.end());
	test(contracted_graph.no_nodes() == 3);

}

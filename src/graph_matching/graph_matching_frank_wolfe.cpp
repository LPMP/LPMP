#include "graph_matching/graph_matching_frank_wolfe.h"
#include "graph_matching/graph_matching_frank_wolfe_impl.hxx"

namespace LPMP {

    graph_matching_frank_wolfe::graph_matching_frank_wolfe(const graph_matching_input& instance, const graph_matching_input::labeling& l, const graph_matching_frank_wolfe_options o)
    {
        if(instance.quadratic_terms.size() >= std::pow(instance.no_left_nodes,2) * std::pow(instance.no_right_nodes,2) / 20) {
            //solver = std::make_unique<graph_matching_frank_wolfe_dense>(instance, l);
            solver = new graph_matching_frank_wolfe_dense(instance, l, o);
        } else {
            //solver = std::make_unique<graph_matching_frank_wolfe_sparse>(instance, l);
            solver = new graph_matching_frank_wolfe_sparse(instance, l, o);
        } 
    }

    graph_matching_input::labeling graph_matching_frank_wolfe::solve()
    {
        return std::visit([](auto&& s){
                s->solve();
                return s->get_solution();
                }, solver);
    }

}

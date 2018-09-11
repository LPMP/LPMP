#include "test_mrf.hxx"
#include "tree_decomposition.hxx"

using namespace LPMP;

int main() 
{
    using FMC = FMC_SRMP;
    using VisitorType = StandardVisitor;

    using SolverType = Solver<LP<FMC>,VisitorType>;
    SolverType solver(solver_options);
    auto& mrf = solver.template GetProblemConstructor<0>();

    // build tree
    auto* u1 = mrf.add_unary_factor(std::vector<REAL>({1.0, 0.0}));
    auto* u2 = mrf.add_unary_factor(std::vector<REAL>({0.0, 2.0}));
    auto* u3 = mrf.add_unary_factor(std::vector<REAL>({0.0, 2.0}));
    auto* u4 = mrf.add_unary_factor(std::vector<REAL>({1.0, 0.0}));

    auto potts_pot = construct_potts(2,2, 0.0, 1.0);
    auto* p12 = mrf.add_pairwise_factor(0,1, potts_pot);
    auto* p23 = mrf.add_pairwise_factor(1,2, potts_pot);
    auto* p34 = mrf.add_pairwise_factor(2,3, potts_pot);

    auto trees = mrf.compute_forest_cover();
    test(trees.size() == 1);
    trees[0].populate_factors();
    test(trees[0].solve() == 2.0);
}

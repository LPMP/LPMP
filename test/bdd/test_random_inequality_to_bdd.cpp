#include "bdd/bdd_min_marginal_averaging.h"
#include "bdd/bdd_anisotropic_diffusion.h"
#include "bdd/convert_pb_to_bdd.h"
#include <vector>
#include <random>
#include "test.h"

// TODO: rename single to random

using namespace LPMP;

// coefficients, inequality type, right_hand_side
std::tuple<std::vector<int>, inequality_type, int> generate_random_inequality(const std::size_t nr_vars)
{
    std::uniform_int_distribution<> d(-10,10);
    std::mt19937 gen;

    std::vector<int> coefficients;
    for(std::size_t i=0; i<nr_vars; ++i)
        coefficients.push_back( d(gen) );

    inequality_type ineq = [&]() {
        const int r = d(gen);
        if(r > 2)
            return inequality_type::smaller_equal;
        else if(r < -2)
            return inequality_type::greater_equal;
        else
            return inequality_type::equal;
    }();

    return {coefficients, ineq, d(gen)}; 
}

std::vector<double> generate_random_costs(const std::size_t nr_vars)
{
    std::uniform_int_distribution<> d(-10,10);
    std::mt19937 gen;

    std::vector<double> coefficients;
    for(std::size_t i=0; i<nr_vars; ++i)
        coefficients.push_back( d(gen) ); 

    return coefficients;
}

template<typename LHS_ITERATOR, typename COST_ITERATOR, typename SOL_ITERATOR>
double min_cost_impl(LHS_ITERATOR lhs_begin, LHS_ITERATOR lhs_end, const inequality_type ineq, const int rhs, COST_ITERATOR cost_begin, COST_ITERATOR cost_end, SOL_ITERATOR sol_begin, const double partial_cost, double& best_current_sol)
{
    assert(std::distance(lhs_begin, lhs_end) == std::distance(cost_begin, cost_end));

    if(lhs_begin == lhs_end) {
        if(ineq == inequality_type::equal) {
            return rhs == 0 ? 0.0 : std::numeric_limits<double>::infinity();
        } else if(ineq == inequality_type::smaller_equal) {
            return rhs >= 0 ? 0.0 : std::numeric_limits<double>::infinity();
        } else if(ineq == inequality_type::greater_equal) {
            return rhs <= 0 ? 0.0 : std::numeric_limits<double>::infinity();
        } 
    }

    const double zero_cost = min_cost_impl(lhs_begin+1, lhs_end, ineq, rhs, cost_begin+1, cost_end, sol_begin+1, partial_cost, best_current_sol);
    const double one_cost = min_cost_impl(lhs_begin+1, lhs_end, ineq, rhs - *lhs_begin, cost_begin+1, cost_end, sol_begin+1, partial_cost + *cost_begin, best_current_sol) + *cost_begin;

    const double sub_tree_cost = std::min(zero_cost, one_cost);
    const double cur_cost = partial_cost + sub_tree_cost;
    if(cur_cost <= best_current_sol) {
        best_current_sol = cur_cost;
        *sol_begin = zero_cost < one_cost ? 0 : 1; 
    }

    return std::min(zero_cost, one_cost);
}

template<typename LHS_ITERATOR, typename COST_ITERATOR>
std::tuple<double, std::vector<char>> min_cost(LHS_ITERATOR lhs_begin, LHS_ITERATOR lhs_end, const inequality_type ineq, const int rhs, COST_ITERATOR cost_begin, COST_ITERATOR cost_end)
{
    std::vector<char> sol(std::distance(lhs_begin, lhs_end));

    double opt_val = std::numeric_limits<double>::infinity();
    const double opt_val_2 = min_cost_impl(lhs_begin, lhs_end, ineq, rhs, cost_begin, cost_end, sol.begin(), 0.0, opt_val);
    assert(opt_val == opt_val_2);

    return {opt_val, sol};
}

template<typename LHS_ITERATOR, typename COST_ITERATOR>
double exp_sum_impl(LHS_ITERATOR lhs_begin, LHS_ITERATOR lhs_end, const inequality_type ineq, const int rhs, COST_ITERATOR cost_begin, COST_ITERATOR cost_end, const double partial_sum)
{
    assert(std::distance(lhs_begin, lhs_end) == std::distance(cost_begin, cost_end));

    if(lhs_begin == lhs_end) {
        if(ineq == inequality_type::equal) {
            return rhs == 0 ? partial_sum : 0.0;
        } else if(ineq == inequality_type::smaller_equal) {
            return rhs >= 0 ? partial_sum : 0.0;
        } else if(ineq == inequality_type::greater_equal) {
            return rhs <= 0 ? partial_sum : 0.0;
        } 
    }

    const double zero_cost = min_cost_impl(lhs_begin+1, lhs_end, ineq, rhs, cost_begin+1, cost_end, sol_begin+1, partial_cost);
    const double one_cost = min_cost_impl(lhs_begin+1, lhs_end, ineq, rhs - *lhs_begin, cost_begin+1, cost_end, sol_begin+1, partial_cost + *cost_begin) + *cost_begin;

    const double sub_tree_sum = zero_cost + one_cost;
    const double cur_cost = partial_cost + sub_tree_cost;

    return std::min(zero_cost, one_cost);
}

template<typename LHS_ITERATOR, typename COST_ITERATOR>
double exp_sum(LHS_ITERATOR lhs_begin, LHS_ITERATOR lhs_end, const inequality_type ineq, const int rhs, COST_ITERATOR cost_begin, COST_ITERATOR cost_end)
{
    const double sum = min_cost_impl(lhs_begin, lhs_end, ineq, rhs, cost_begin, cost_end, sol.begin(), 0.0);
    return sum;
} 

template<typename BDD_SOLVER>
void test_random_inequality()
{
    Cudd bdd_mgr;
    bdd_converter converter(bdd_mgr);

    for(std::size_t nr_vars = 3; nr_vars <= 15; ++nr_vars) {
        const auto [coefficients, ineq, rhs] = generate_random_inequality(nr_vars);
        for(const auto c : coefficients) {
            std::cout << c << " ";
        }
        if(ineq == inequality_type::equal)
            std::cout << " = ";
        if(ineq == inequality_type::smaller_equal)
            std::cout << " <= ";
        if(ineq == inequality_type::greater_equal)
            std::cout << " >= ";
        std::cout << rhs << "\n";

        auto bdd = converter.convert_to_bdd(coefficients.begin(), coefficients.end(), ineq, rhs);
        if(bdd.nodeCount() < 2) 
            continue;
        BDD_SOLVER bdds;
        std::vector<std::size_t> vars(nr_vars);
        std::iota (std::begin(vars), std::end(vars), 0);
        bdds.add_bdd(bdd, vars.begin(), vars.end(), bdd_mgr);
        //bdds.export_dot(std::cout);
        bdds.init(); 
        const std::vector<double> costs = generate_random_costs(nr_vars);
        std::cout << "cost: ";
        for(const auto x : costs)
            std::cout << x << " ";
        std::cout << "\n"; 
        bdds.set_costs(costs.begin(), costs.end());
        const double backward_lb = bdds.lower_bound();
        const auto [enumeration_lb, sol] = min_cost(coefficients.begin(), coefficients.end(), ineq, rhs, costs.begin(), costs.end());
        std::cout << "enumeration lb = " << enumeration_lb << ", backward lb = " << backward_lb << "\n";
        test(std::abs(backward_lb - enumeration_lb) <= 1e-8);
        std::cout << "cost of primal = " << bdds.evaluate(sol.begin(), sol.end()) << "\n";
        std::cout << "primal size = " << sol.size() << "\n";
        for(const auto x : sol)
            std::cout << int(x) << " ";
        std::cout << "\n";
        test(std::abs(enumeration_lb - bdds.evaluate(sol.begin(), sol.end())) <= 1e-8);
    } 
}

int main(int argc, char** arv)
{
    test_random_inequality<bdd_min_marginal_averaging>();
    test_random_inequality<bdd_anisotropic_diffusion>();
}

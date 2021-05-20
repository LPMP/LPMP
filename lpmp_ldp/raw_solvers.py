import numpy as np
import bindings.graph_matching_py as gm
import bindings.multigraph_matching_py as mgm


def gm_solver(costs, quadratic_costs, edges_left, edges_right, solver_params, verbose=False):
    """
    A thin python wrapper of the solver Graph Matching solver. Computes min-cost matching between two directed graphs
    G_1 = (V_1, E_1) and G_2 = (V_2, E_2).

    @param costs: np.array of shape (|V_1|, |V_2|) with unary matching costs
    @param quadratic_costs: np.array of shape (|E_1|, |E_2|) with pairwise matching costs
    @param edges_left: np.array of shape (|E_1|, 2) with edges of G_1 (pairs of vertex indices)
    @param edges_right: np.array of shape (|E_2|, 2) with edges of G_2 (pairs of vertex indices)
    @param solver_params: dict of command line flags to pass to the solver (see solver documentation)
    @param verbose: bool, if true print summary of input data AND raw solver output
    @return: np.array of shape (|V_1|, |V_2|) with 0/1 values capturing the min-cost matching produced by the solver,
             np.array of shape (|E_1|, |E_2|) with 0/1 values capturing which pairwise costs were paid in suggested
             matching
    """
    if verbose:
        print(
            f"Building graph with\t"
            f"#vertices left: {costs.shape[0]}\t"
            f"#vertices right: {costs.shape[1]}\t"
            f"#edges left: {len(edges_left)}\t"
            f"#edges right: {len(edges_right)}\t"
            f"#quadratic terms: {quadratic_costs.size}"
        )

    instance = gm.graph_matching_input(
        costs, quadratic_costs, np.array(edges_left, dtype=np.intc), np.array(edges_right, dtype=np.intc)
    )

    params = ["tmp", f"-v {int(verbose)}"] + [f"--{key} {val}" for key, val in solver_params.items()]
    solver = gm.graph_matching_message_passing_solver(params)

    solver.construct(instance)
    solver.solve()

    result = solver.result()
    costs_paid, quadratic_costs_paid = result.result_masks(
        costs, quadratic_costs, np.array(edges_left), np.array(edges_right)
    )

    return costs_paid, quadratic_costs_paid


def mgm_solver(unary_costs, quadratic_costs, edges, solver_params, verbose=False):
    """
    A thin python wrapper of the solver Multigraph Matching solver. Computes min-cost matching of k directed graphs
    G_1 = (V_1, E_1),...,G_k = (V_k, E_k). The suggested matching is "globally consistent" i.e. obey cycle consistency

    @param unary_costs: a list of {k choose 2} unary costs - np.arrays of shape (|V_i|, |V_j|)
    for i \neq j in lexicographical order
    @param quadratic_costs: a list of {k choose 2} pairwise costs - np.arrays of shape (|E_i|, |E_j|) for i \neq j in
    lexicographical order
    @param edges: a list of k edge descriptions - np.arrays of shape (|E_i|, 2)
    @param solver_params: dict of command line flags to pass to the solver (see solver documentation)
    @param verbose: bool, if true print raw solver output
    @return: list of {k choose 2} np.arrays of shape (|V_i|, |V_j|) with 0/1 values capturing the min-cost matching
             list of {k choose 2} np.array of shape (|E_i|, |E_j|) with 0/1 values capturing which pairwise costs were
             paid in suggested matching
    """
    instance = mgm.multigraph_matching_input(unary_costs, quadratic_costs, edges)

    params = ["", f"-v {int(verbose)}"] + [
        f"--{key} {val}" if val else f"--{key}" for key, val in solver_params.items()
    ]
    solver = mgm.multigraph_matching_message_passing_solver(params)

    solver.construct(instance)
    solver.solve()
    result = solver.result()

    costs_paid, quadratic_costs_paid = result.result_masks(unary_costs, quadratic_costs, edges)
    return costs_paid, quadratic_costs_paid


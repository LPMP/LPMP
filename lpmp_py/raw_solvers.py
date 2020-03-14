import numpy as np
import bindings.graph_matching_py as gm
import bindings.multigraph_matching_py as mgm


def gm_solver(costs, quadratic_costs, edges_left, edges_right, solver_params, verbose=False):
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


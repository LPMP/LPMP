import numpy as np
import torch

from utils.utils import torch_to_numpy_list, numpy_to_torch_list
from utils.utils import lexico_iter_pairs
import bindings.multigraph_matching_py as mgm


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


def detach(list_of_tensors):
    return [x.detach() for x in list_of_tensors]


class MultiGraphMatchingSolver(torch.autograd.Function):
    @staticmethod
    def forward(ctx, params, *all_costs_tup):
        l = len(all_costs_tup) // 2
        costs_l, quadratic_costs_l = all_costs_tup[:l], all_costs_tup[l:]
        assert len(costs_l) == len(quadratic_costs_l)
        costs_paid_l, quadratic_costs_paid_l = mgm_solver(
            unary_costs=torch_to_numpy_list(costs_l),
            quadratic_costs=torch_to_numpy_list(quadratic_costs_l),
            edges=torch_to_numpy_list(params["edges"]),
            solver_params=params["solver_params"],
        )
        device = costs_l[0].device
        costs_paid_l = numpy_to_torch_list(costs_paid_l, device=device, dtype=torch.float32)
        quadratic_costs_paid_l = numpy_to_torch_list(quadratic_costs_paid_l, device=device, dtype=torch.float32)

        ctx.params = params
        ctx.costs_l = detach(costs_l)
        ctx.costs_paid_l = detach(costs_paid_l)
        ctx.quadratic_costs_l = detach(quadratic_costs_l)
        ctx.quadratic_costs_paid_l = detach(quadratic_costs_paid_l)
        return (*costs_paid_l, *quadratic_costs_paid_l)

    @staticmethod
    def backward(ctx, *all_grad_costs_paid_tup):
        l = len(all_grad_costs_paid_tup) // 2
        grad_costs_paid_l, grad_quadratic_costs_paid_l = all_grad_costs_paid_tup[:l], all_grad_costs_paid_tup[l:]
        params = ctx.params
        lambda_val = params["lambda_val"]
        epsilon_val = 1e-8

        all_costs_prime = [
            costs + lambda_val * grad_costs_paid for costs, grad_costs_paid in zip(ctx.costs_l, grad_costs_paid_l)
        ]
        all_quadratic_costs_prime = [
            quadratic_costs + lambda_val * grad_quadratic_costs_paid
            for quadratic_costs, grad_quadratic_costs_paid in zip(ctx.quadratic_costs_l, grad_quadratic_costs_paid_l)
        ]

        all_costs_paid_prime, all_quadratic_costs_paid_prime = mgm_solver(
            unary_costs=torch_to_numpy_list(all_costs_prime),
            quadratic_costs=torch_to_numpy_list(all_quadratic_costs_prime),
            edges=torch_to_numpy_list(params["edges"]),
            solver_params=params["solver_params"],
        )
        device = grad_costs_paid_l[0].device
        all_costs_paid_prime = numpy_to_torch_list(all_costs_paid_prime, device=device, dtype=torch.float32)
        all_quadratic_costs_paid_prime = numpy_to_torch_list(
            all_quadratic_costs_paid_prime, device=device, dtype=torch.float32
        )

        grad_costs_l = [
            -(costs_paid - costs_paid_prime) / (lambda_val + epsilon_val)
            for costs_paid, costs_paid_prime in zip(ctx.costs_paid_l, all_costs_paid_prime)
        ]
        grad_quadratic_costs_l = [
            -(quadratic_costs_paid - quadratic_costs_paid_prime) / (lambda_val + epsilon_val)
            for quadratic_costs_paid, quadratic_costs_paid_prime in zip(
                ctx.quadratic_costs_paid_l, all_quadratic_costs_paid_prime
            )
        ]

        return (None, *grad_costs_l, *grad_quadratic_costs_l)


class MultiGraphMatchingModule(torch.nn.Module):
    def __init__(self, edges_batch_list, num_vertices_batch_list, lambda_val, solver_params):
        super().__init__()
        self.solver = MultiGraphMatchingSolver()
        self.edges_batch_list = edges_batch_list
        self.num_vertices_batch_list = num_vertices_batch_list
        self.params = {"lambda_val": lambda_val, "solver_params": solver_params}

    def forward(self, costs_batch_list, quadratic_costs_batch_list):
        def params_generator():
            for edges in zip(*self.edges_batch_list):
                yield {"edges": [e.T for e in edges], **self.params}

        def costs_generator():
            zipped = zip(
                zip(*self.edges_batch_list),
                zip(*self.num_vertices_batch_list),
                zip(*costs_batch_list),
                zip(*quadratic_costs_batch_list),
            )

            for edges_l, num_vertices_l, costs_l, quadratic_costs_l in zipped:
                truncated_costs_l = []
                truncated_quadratic_costs_l = []
                leftover_costs_l = []
                leftover_quadratic_costs_l = []

                for (
                    costs,
                    quadratic_costs,
                    (edges_s, edges_t),
                    (num_vertices_s, num_vertices_t),
                ) in zip(costs_l, quadratic_costs_l, lexico_iter_pairs(edges_l), lexico_iter_pairs(num_vertices_l)):
                    truncated_costs = costs[:num_vertices_s, :num_vertices_t]
                    assert quadratic_costs.shape[0] == edges_s.shape[1], (quadratic_costs.shape, edges_s.shape)
                    assert quadratic_costs.shape[1] == edges_t.shape[1], (quadratic_costs.shape, edges_t.shape)
                    truncated_quadratic_costs = quadratic_costs[: edges_s.shape[-1], : edges_t.shape[-1]]
                    leftover_costs = (truncated_costs.abs().sum() - costs.abs().sum()).abs()
                    assert leftover_costs < 1e-5, leftover_costs
                    leftover_quadratic_costs = (
                        truncated_quadratic_costs.abs().sum() - quadratic_costs.abs().sum()
                    ).abs()
                    assert leftover_quadratic_costs < 1e-5, leftover_quadratic_costs

                    truncated_costs_l.append(truncated_costs)
                    truncated_quadratic_costs_l.append(truncated_quadratic_costs)
                    leftover_costs_l.append(leftover_costs)
                    leftover_quadratic_costs_l.append(leftover_quadratic_costs)
                yield truncated_costs_l, truncated_quadratic_costs_l

        batch_size = len(costs_batch_list[0])
        result = []
        for costs_batch in costs_batch_list:
            max_dimension_x = max(x.shape[0] for x in costs_batch)
            max_dimension_y = max(x.shape[1] for x in costs_batch)
            result.append(torch.zeros(size=(batch_size, max_dimension_x, max_dimension_y)).cuda())

        for i, (params_l, (un_costs_l, quad_costs_l), num_vertices_l) in enumerate(
            zip(params_generator(), costs_generator(), zip(*self.num_vertices_batch_list))
        ):
            all_costs_paid_tup = self.solver.apply(params_l, *(un_costs_l + quad_costs_l))
            un_costs_paid_tup = all_costs_paid_tup[: len(all_costs_paid_tup) // 2]
            for j, ((num_vertices_s, num_vertices_t), costs_paid) in enumerate(
                zip(lexico_iter_pairs(num_vertices_l), un_costs_paid_tup)
            ):
                result[j][i, :num_vertices_s, :num_vertices_t] = costs_paid  # Only unary matching returned
        return result

import torch

from .utils import torch_to_numpy_list, numpy_to_torch_list, lexico_iter_pairs, detach
from ..raw_solvers import mgm_solver

class MultiGraphMatchingSolver(torch.autograd.Function):
    """
        Multigraph Matching solver as a torch.Function where the backward pass is provided by
        [1] 'Vlastelica* M, Paulus* A., Differentiation of Blackbox Combinatorial Solvers, ICLR 2020'
    """
    @staticmethod
    def forward(ctx, params, *all_costs_tup):
        """
            Implementation of the forward pass of min-cost matching between k directed graphs
            G_1 = (V_1, E_1),..., G_k = (V_k, E_k)

            @param ctx: context for backpropagation
            @param params: a dict of additional params. Must contain:
                    edges: a list of k torch.tensor with shapes (|E_i|, 2) describing edges of G_i,
                    lambda_val: float/np.float32/torch.float32, the value of lambda for computing the gradient with [1]
                    solver_params: a dict of command line parameters to the solver (see solver documentation)
            @param all_costs_tup: a list of 2*{k choose 2} torch.Tensors:
                    The first {k choose 2} describe unary costs - torch.Tensors of shape (|V_i|, |V_j|) for i \neq j
                        in lexicographical order.
                    The second {k choose 2} describe quadratic costs - torch.Tensors of shape (|E_i|, |E_j|) for
                        i \neq j in lexicographical order
            @return: a tuple of list of 2*{k choose 2} torch.Tensors:
                    The first {k choose 2} are paid unary costs - torch.Tensors of shape (|V_i|, |V_j|) with 0/1 values
                    capturing the min-cost matching.
                    The second {k choose 2} describe paid quadratic costs torch.Tensors of shape (|E_i|, |E_j|)
                    with 0/1 values capturing which pairwise costs were paid in suggested matching.
        """
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
    """
        Torch module for handling batches of Multigraph Matching Instances.
    """
    def __init__(self, edges_batch_list, num_vertices_batch_list, lambda_val, solver_params):
        """
        Prepares a module for a batch of k multigraph matching instances each with n graphs.

        @param edges_batch_list: A list of k lists, each containing n torch.Tensors of shape (num_edges, 2) describing
        edges of all graphs in all instances.
        @param num_vertices_batch_list: A list of k lists, each containing n integers, the numbers of vertices in the
        participating graphs
        @param lambda_val: lambda value for backpropagation by [1]
        @param solver_params: a dict of command line parameters to the solver (see solver documentation)
        """
        super().__init__()
        self.solver = MultiGraphMatchingSolver()
        self.edges_batch_list = edges_batch_list
        self.num_vertices_batch_list = num_vertices_batch_list
        self.params = {"lambda_val": lambda_val, "solver_params": solver_params}

    def forward(self, costs_batch_list, quadratic_costs_batch_list):
        """
        Forward pass for a batch of k multigraph matching instances on n graphs
        @param costs_batch_list: a list of k lists, each with {n choose 2} unary cost torch.Tensors
        @param quadratic_costs_batch_list: a list of k lists, each with {n choose 2} pariwise cost torch.Tensors
        @return: a list of k lists, each with {n choose 2} unary cost torch.Tensors with the suggested matchings
        """
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

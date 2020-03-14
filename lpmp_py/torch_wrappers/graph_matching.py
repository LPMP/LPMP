import torch
from ..raw_solvers import gm_solver


class GraphMatchingSolver(torch.autograd.Function):
    @staticmethod
    def forward(ctx, costs, quadratic_costs, params):
        device = costs.device
        costs_paid, quadratic_costs_paid = gm_solver(
            costs=costs.cpu().detach().numpy(),
            quadratic_costs=quadratic_costs.cpu().detach().numpy(),
            edges_left=params["edges_left"].cpu().detach().numpy(),
            edges_right=params["edges_right"].cpu().detach().numpy(),
            solver_params=params["solver_params"],
        )
        costs_paid = torch.from_numpy(costs_paid).to(torch.float32).to(device)
        quadratic_costs_paid = torch.from_numpy(quadratic_costs_paid).to(torch.float32).to(device)
        ctx.params = params
        ctx.save_for_backward(costs, costs_paid, quadratic_costs, quadratic_costs_paid)
        return costs_paid, quadratic_costs_paid

    @staticmethod
    def backward(ctx, grad_costs_paid, grad_quadratic_costs_paid):
        costs, costs_paid, quadratic_costs, quadratic_costs_paid = ctx.saved_tensors
        device = costs.device
        lambda_val = ctx.params["lambda_val"]
        epsilon_val = 1e-8
        assert grad_costs_paid.shape == costs.shape and grad_quadratic_costs_paid.shape == quadratic_costs.shape

        costs_prime = costs + lambda_val * grad_costs_paid
        quadratic_costs_prime = quadratic_costs + lambda_val * grad_quadratic_costs_paid
        costs_paid_prime, quadratic_costs_paid_prime = gm_solver(
            costs=costs_prime.cpu().detach().numpy(),
            quadratic_costs=quadratic_costs_prime.cpu().detach().numpy(),
            edges_left=ctx.params["edges_left"].cpu().detach().numpy(),
            edges_right=ctx.params["edges_right"].cpu().detach().numpy(),
            solver_params=ctx.params["solver_params"],
        )
        costs_paid_prime = torch.from_numpy(costs_paid_prime).to(torch.float32).to(device)
        quadratic_costs_paid_prime = torch.from_numpy(quadratic_costs_paid_prime).to(torch.float32).to(device)

        grad_costs = -(costs_paid - costs_paid_prime) / (lambda_val + epsilon_val)
        grad_quadratic_costs = -(quadratic_costs_paid - quadratic_costs_paid_prime) / (lambda_val + epsilon_val)

        return grad_costs, grad_quadratic_costs, None


class GraphMatchingModule(torch.nn.Module):
    def __init__(
        self,
        edges_left_batch,
        edges_right_batch,
        num_vertices_s_batch,
        num_vertices_t_batch,
        lambda_val,
        solver_params,
    ):
        super().__init__()
        self.solver = GraphMatchingSolver()
        self.edges_left_batch = edges_left_batch
        self.edges_right_batch = edges_right_batch
        self.num_vertices_s_batch = num_vertices_s_batch
        self.num_vertices_t_batch = num_vertices_t_batch
        self.params = {"lambda_val": lambda_val, "solver_params": solver_params}

    def forward(self, costs_batch, quadratic_costs_batch):
        def params_generator():
            for edges_left, edges_right in zip(self.edges_left_batch, self.edges_right_batch):
                yield {"edges_left": edges_left.T, "edges_right": edges_right.T, **self.params}

        def costs_generator():
            zipped = zip(
                self.edges_left_batch,
                self.edges_right_batch,
                self.num_vertices_s_batch,
                self.num_vertices_t_batch,
                costs_batch,
                quadratic_costs_batch,
            )
            for edges_left, edges_right, num_vertices_s, num_vertices_t, costs, quadratic_costs in zipped:
                truncated_costs = costs[:num_vertices_s, :num_vertices_t]
                assert quadratic_costs.shape[0] == edges_left.shape[1], (quadratic_costs.shape, edges_left.shape)
                assert quadratic_costs.shape[1] == edges_right.shape[1], (quadratic_costs.shape, edges_right.shape)
                truncated_quadratic_costs = quadratic_costs[: edges_left.shape[-1], : edges_right.shape[-1]]
                leftover_costs = (truncated_costs.abs().sum() - costs.abs().sum()).abs()
                assert leftover_costs < 1e-5, leftover_costs
                leftover_quadratic_costs = (truncated_quadratic_costs.abs().sum() - quadratic_costs.abs().sum()).abs()
                assert leftover_quadratic_costs < 1e-5, leftover_quadratic_costs
                yield truncated_costs, truncated_quadratic_costs

        batch_size = len(costs_batch)
        max_dimension_x = max(x.shape[0] for x in costs_batch)
        max_dimension_y = max(x.shape[1] for x in costs_batch)
        result = torch.zeros(size=(batch_size, max_dimension_x, max_dimension_y)).cuda()
        for i, (params, costs, num_vertices_s, num_vertices_t) in enumerate(
            zip(params_generator(), costs_generator(), self.num_vertices_s_batch, self.num_vertices_t_batch)
        ):
            tmp = self.solver.apply(costs[0], costs[1], params)
            result[i, :num_vertices_s, :num_vertices_t] = tmp[0]  # Only unary matching returned
        return result

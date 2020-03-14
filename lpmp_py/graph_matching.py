import bindings.graph_matching_py as gm
import numpy as np
import torch


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


class GraphMatching(torch.nn.Module):
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


# testing


def random_instance(size_left, size_right):
    graph_left = np.ones([size_left] * 2, dtype=int)
    graph_left = graph_left - np.diag(graph_left.diagonal())

    it = np.nditer(graph_left, flags=["multi_index"])
    edges_left = [it.multi_index for i in it if i == 1]

    graph_right = np.random.uniform(0, 1, size=(size_right, size_right)) > 0.5
    graph_right = np.triu(graph_right, 1)

    it = np.nditer(graph_right, flags=["multi_index"])
    edges_right = [it.multi_index for i in it if i == 1]

    costs_unary = 0.1 * np.array(2.0 * np.random.randint(-10, -1, [size_left, size_right]) - 1.0, dtype=np.double)

    costs_binary = 0.1 * np.array(
        2.0 * np.random.randint(-15, -1, [len(edges_left), len(edges_right)]) - 1.0, dtype=np.double
    )

    return edges_left, edges_right, costs_unary, costs_binary


def random_instance_batch(max_num_vertices, batch_size):
    num_vertices_batch = np.random.randint(4, max_num_vertices, [batch_size])

    instances = [random_instance(num_vertices, num_vertices) for num_vertices in num_vertices_batch]

    edges_left_batch, edges_right_batch, costs_unary_batch, costs_binary_batch = zip(*instances)

    max_num_vertices = max(num_vertices_batch)
    max_edges_left = max([len(edges_left) for edges_left in edges_left_batch])
    max_edges_right = max([len(edges_right) for edges_right in edges_right_batch])

    costs_unary_batch_numpy = np.zeros([batch_size, max_num_vertices, max_num_vertices])
    for i, (unary_costs, num_vertices) in enumerate(zip(costs_unary_batch, num_vertices_batch)):
        costs_unary_batch_numpy[i, :num_vertices, :num_vertices] = unary_costs

    costs_binary_batch_numpy = np.zeros([batch_size, max_edges_left, max_edges_right])
    for i, (binary_costs, edges_left, edges_right) in enumerate(
        zip(costs_binary_batch, edges_left_batch, edges_right_batch)
    ):
        costs_binary_batch_numpy[i, : len(edges_left), : len(edges_right)] = binary_costs

    return (
        (edges_left_batch, edges_right_batch, num_vertices_batch),
        (costs_unary_batch_numpy, costs_binary_batch_numpy),
    )


if __name__ == "__main__":

    verbose = True

    graph_info_batch, costs_batch = random_instance_batch(max_num_vertices=15, batch_size=8)

    solver_params = {"timeout": 200, "primalComputationInterval": 1, "maxIter": 30}

    Solver = GraphMatching(*graph_info_batch, 1.0, solver_params)

    costs_unary, costs_binary = costs_batch
    costs_unary_torch = torch.from_numpy(costs_unary)
    costs_unary_torch.requires_grad = True
    costs_binary_torch = torch.from_numpy(costs_binary)
    costs_binary_torch.requires_grad = True

    import time

    beg = time.time()
    out_torch = Solver(costs_unary_torch, costs_binary_torch)

    random_mask = torch.randn_like(out_torch)
    faked_loss = (out_torch * random_mask).sum()
    faked_loss.backward()

    secs = time.time() - beg
    print("Unary gradient:", costs_unary_torch.grad.abs().sum())
    print("Quadratic gradient: ", costs_binary_torch.grad.abs().sum())
    print(f"Time elapsed: {secs}")

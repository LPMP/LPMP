import torch
import torch.multiprocessing as mp
import numpy as np
from .raw_solvers import mc_solver_gaec, mwc_solver
from .utils import run_func

def solve_mc(batch_index, return_dict, instance):
    node_labels, edge_labels = mc_solver_gaec(instance['edge_indices'], instance['edge_costs'])
    return_dict[batch_index] = {'node_labels': node_labels, 'edge_labels': edge_labels}

def solve_mwc(batch_index, return_dict, instance):
    node_labels, edge_labels, sol_cost = mwc_solver(instance['edge_indices'], instance['edge_costs'], instance['node_costs'])
    return_dict[batch_index] = {'node_labels': node_labels, 'edge_labels': edge_labels}

def convert_contiguous(cluster_labels):
    unique_labels = np.unique(cluster_labels)
    mapping = np.zeros((unique_labels.max() + 1, 1))
    mapping[unique_labels] = 1
    mapping = (np.cumsum(mapping) - 1).astype(np.int32)
    def contiguous_mapping_func(entry):
        return mapping[entry]
    contiguous_mapping = np.vectorize(contiguous_mapping_func)
    cluster_labels_cont = contiguous_mapping(cluster_labels)
    return cluster_labels_cont

def map_edge_indices(edge_indices, node_mapping):
    row = edge_indices[:, 0]
    col = edge_indices[:, 1]
    def node_mapping_func(entry):
        return node_mapping[entry]
    node_mapping_vect = np.vectorize(node_mapping_func)
    row = node_mapping_vect(row)
    col = node_mapping_vect(col)
    diff_cc_mask = row != col
    return np.stack((row[diff_cc_mask], col[diff_cc_mask]), 1)

def get_complete_graph(num_nodes, device):
    row, col = torch.meshgrid(torch.arange(num_nodes), torch.arange(num_nodes), device = device)
    row = row.flatten()
    col = col.flatten()
    non_diag = row != col
    return torch.stack((row[non_diag], col[non_diag]), 1)
    
class MulticutSolver(torch.autograd.Function):
    @staticmethod
    def run_multicut_batch(edge_indices, edge_costs):
        device = edge_costs.device
        batch_size = edge_costs.shape[0]
        assert batch_size == edge_indices.shape[0]
        edge_indices_np = edge_indices.detach().cpu().numpy()
        edge_costs_np = edge_costs.detach().cpu().numpy()
        batched_instances = []
        for b in range(batch_size):
            batched_instances.append({'edge_indices': edge_indices_np[b, ::], 'edge_costs': edge_costs_np[b, ::]})

        results = run_func(solve_mc, batched_instances)
        node_labels_batch = []
        edge_labels_batch = []
        new_edge_indices_batch = []
        for b in sorted(results.keys()):
            node_labels_cont = convert_contiguous(results[b]['node_labels'])
            new_edge_indices_batch.append(torch.from_numpy(map_edge_indices(edge_indices_np[b, ::], node_labels_cont)).to(device).long())
            # new_edge_indices_batch.append(get_complete_graph(num_nodes, device))
            node_labels_one_hot = torch.nn.functional.one_hot(torch.from_numpy(node_labels_cont).to(device).long())
            node_labels_batch.append(node_labels_one_hot.to(torch.float32).to(device))
            
            edge_labels_batch.append(results[b]['edge_labels'])

        edge_labels_batch = torch.from_numpy(np.stack(edge_labels_batch, 0)).to(torch.float32).to(device)

        return node_labels_batch, edge_labels_batch, new_edge_indices_batch

    @staticmethod
    def run_multiway_cut_batch(edge_indices, edge_costs, node_costs):
        device = edge_costs.device
        batch_size = edge_costs.shape[0]
        assert batch_size == edge_indices.shape[0]
        assert batch_size == node_costs.shape[0]
        edge_indices_np = edge_indices.detach().cpu().numpy()
        node_costs_np = node_costs.detach().cpu().numpy()
        edge_costs_np = edge_costs.detach().cpu().numpy()
        batched_instances = []
        for b in range(batch_size):
            batched_instances.append({'node_costs': node_costs_np[b, ::], 'edge_indices': edge_indices_np[b, ::], 'edge_costs': edge_costs_np[b, ::]})

        results = run_func(solve_mwc, batched_instances)
        edge_labels_batch = []
        for b in sorted(results.keys()):
            edge_labels_batch.append(results[b]['edge_labels'])
        edge_labels_batch = torch.from_numpy(np.stack(edge_labels_batch, 0)).to(torch.float32).to(device)
        return edge_labels_batch

    @staticmethod
    def run_multiway_cut_perturb_batch(edge_indices, edge_costs, grad_node_labels, orig_edge_labels, lambda_val_start, lambda_val_end, num_grad_samples):
        device = edge_costs.device
        batch_size = edge_costs.shape[0]
        assert batch_size == edge_indices.shape[0]
        assert batch_size == len(grad_node_labels)
        sampled_lambdas = torch.rand((num_grad_samples), device=device) * (lambda_val_end - lambda_val_start) + lambda_val_start
        
        edge_indices_np = edge_indices.detach().cpu().numpy()
        edge_costs_np = edge_costs.detach().cpu().numpy()
        batched_instances = []
        for b in range(batch_size):
            for s in range(num_grad_samples):
                current_pert_node_costs_forward = (grad_node_labels[b] * sampled_lambdas[s]).detach().cpu().numpy()
                batched_instances.append({'node_costs': current_pert_node_costs_forward, 'edge_indices': edge_indices_np[b, ::], 'edge_costs': edge_costs_np[b, ::]})

        results = run_func(solve_mwc, batched_instances)
        index = 0
        grad_edge_costs = torch.zeros_like(edge_costs)
        for b in range(batch_size):
            for s in range(num_grad_samples):
                forward_grad = torch.from_numpy(results[index]['edge_labels']).to(torch.float32).to(device)
                index = index + 1
                grad_edge_costs[b, ::] += (forward_grad - orig_edge_labels[b]) / sampled_lambdas[s]

        grad_edge_costs /= num_grad_samples
        return grad_edge_costs

    @staticmethod
    def forward(ctx, edge_costs, params):
        """
        Implementation of the forward pass of multicut on multiple undirected graphs {G_1 = (V_1, E_1), ..., G_B = (V_B, E_B)}
        @param ctx: context for backpropagation
        @param edge_costs: torch.Tensor of shape (|B| x |E|) containing edge costs.
        @param params: a dict of additional params. Must contain:
                edge_indices: torch.Tensor of shape (|B| x |E| x 2) containing edge indices.
                lambda_val_start: float, Minimum pertubation
                lambda_val_end: float, Maximum pertubation
                num_grad_samples: Number of gradient samples to average. In each sample a random
                    lambda is sampled in [lambda_val_start, lambda_val_end]
                grad_clusters: Bool. False (Default): Differentiate w.r.t edge labels.
                                     True: Differentiate w.r.t cluster labels. 
        @return: 
            tuples of size |B| where each entry is clustering matrix of shape (|V| x |K_i|), where |K_i| can be unequal to |K_j|.
            torch.Tensor of shape (|B| x |E|) with 0/1 values of edge labeling. 1 means edge is cut.
            torch.Tensor of shape (|B| x |E'| x 2) containing new edge indices after clustering the nodes.
        """
        node_labels_batch, edge_labels_batch, new_edge_indices = MulticutSolver.run_multicut_batch(params['edge_indices'], edge_costs)
        node_labels_batch = tuple(node_labels_batch)
        [ctx.mark_non_differentiable(e) for e in new_edge_indices]
        new_edge_indices = tuple(new_edge_indices)
        if params.get('grad_clusters', False):
            ctx.mark_non_differentiable(edge_labels_batch)
        else:
            ctx.mark_non_differentiable(node_labels_batch)
        ctx.save_for_backward(edge_costs, edge_labels_batch)
        ctx.params = params
        ctx.batch_size = edge_costs.shape[0]

        return node_labels_batch + (edge_labels_batch, ) + new_edge_indices

    @staticmethod
    def backward(*grads):
        """
        Backward pass computation.

        @param ctx: context from the forward pass
        @param grad_node_labels: "dL / d node_labels_batch" tuples of size |B| where each entry is clustering matrix of shape (|V| x |K_i|).
            Where |K_i| can be unequal to |K_j|.
        @param grad_edge_labels: "dL / d edge_labels_batch" torch.Tensor of shape (|B| x |E|).
        @params grad_edge_indices: Non-differentiable.

        @return: dL / edge_costs
        """
        ctx = grads[0]
        batch_size = ctx.batch_size
        params = ctx.params
        edge_costs, edge_labels_batch = ctx.saved_tensors
        device = edge_costs.device
        assert(params.get('grad_clusters', False))
        grad_node_labels = grads[1: 1 + batch_size]
        grad_edge_costs = MulticutSolver.run_multiway_cut_perturb_batch(params['edge_indices'], edge_costs, 
                                                                    grad_node_labels, edge_labels_batch,
                                                                    params['lambda_val_start'], params['lambda_val_end'], params['num_grad_samples'])

        return grad_edge_costs, None

class MulticutModule(torch.nn.Module):
    """
    Torch module for handling batches of Multicut Instances
    """
    def __init__(self, lambda_val_start, lambda_val_end, num_grad_samples, grad_clusters = False):
        """
        @param lambda_val_start: float, Minimum pertubation
        @param lambda_val_end: float, Maximum pertubation
        @param num_grad_samples: Number of gradient samples to average. In each sample a random
            pertubation is sampled in [lambda_val_start, lambda_val_end]
        @param grad_clusters: True: Use cluster labels for backward pass.
        """
        super().__init__()
        self.solver = MulticutSolver()
        self.params = {"lambda_val_start": lambda_val_start, 
                        "lambda_val_end": lambda_val_end, 
                        "num_grad_samples": num_grad_samples, 
                        "grad_clusters": grad_clusters}

    def forward(self, edge_costs, edge_indices):
        self.params['edge_indices'] = edge_indices
        return self.solver.apply(edge_costs, self.params)
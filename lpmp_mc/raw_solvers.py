import numpy as np
import bindings.multicut_py as mc
import bindings.asymmetric_multiway_cut_py as amwc
import bindings.multiway_cut_py as mwc

def mc_solver_gaec(edge_indices, edge_costs):
    mc_instance = mc.multicut_instance(edge_indices, edge_costs)
    edge_labels, node_labels = mc.multicut_gaec(mc_instance)
    edge_labels, node_labels = mc_instance.result_mask(edge_labels, node_labels)
    return node_labels, edge_labels

def mc_solver_gef(edge_indices, edge_costs):
    mc_instance = mc.multicut_instance(edge_indices, edge_costs)
    edge_labels, node_labels = mc.multicut_gef(mc_instance)
    edge_labels, node_labels = mc_instance.result_mask(edge_labels, node_labels)
    return node_labels, edge_labels

def mc_solver_bec(edge_indices, edge_costs):
    mc_instance = mc.multicut_instance(edge_indices, edge_costs)
    edge_labels, node_labels = mc.multicut_bec(mc_instance)
    edge_labels, node_labels = mc_instance.result_mask(edge_labels, node_labels)
    return node_labels, edge_labels

def amwc_solver(node_costs, edge_indices, edge_costs, partitionable = None):
    if partitionable is None:
        partitionable = np.array([], dtype='bool') 
    amc_instance = amwc.asymmetric_multiway_cut_instance(edge_indices, edge_costs, node_costs, partitionable)
    amc_sol = amwc.asymmetric_multiway_cut_gaec(amc_instance)
    edge_labels, node_labels, node_instance_ids = amc_instance.result_mask(amc_sol)
    solver_cost = amc_instance.evaluate(amc_sol)
    return node_labels.astype(np.uint8), node_instance_ids.astype(np.int32), edge_labels.astype(np.uint8), solver_cost

def mwc_solver(edge_indices, edge_costs, node_costs):
    mc_instance = mwc.multiway_cut_instance(edge_indices, edge_costs, node_costs)
    amc_sol = mwc.multiway_cut_gaec(mc_instance)
    edge_labels, node_labels = mc_instance.result_mask(amc_sol)
    solver_cost = mc_instance.evaluate(amc_sol)
    return node_labels.astype(np.uint8), edge_labels.astype(np.uint8), solver_cost
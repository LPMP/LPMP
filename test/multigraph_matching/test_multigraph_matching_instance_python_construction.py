import multigraph_matching_py as mgm
import numpy as np
import random


def construct_mgm_problem(nr_nodes):
    nr_graphs = len(nr_nodes)
    nr_problems = (nr_graphs*(nr_graphs-1))/2

    edges = []
    for n in range(nr_graphs):
        nr_edges = int(nr_nodes[n]*(nr_nodes[n]-1)/2)
        e = np.zeros([nr_edges,2],dtype=np.int)
        c = 0
        for i in range(int(nr_nodes[n])):
            for j in range(i+1,int(nr_nodes[n])):
                e[c,0] = i;
                e[c,1] = j;
                c = c+1 
        edges.append(e)

    assignments = []
    quadratic_terms = []
    for i in range(nr_graphs):
        for j in range(i+1,nr_graphs):
            a = np.random.rand(nr_nodes[i], nr_nodes[j]) - 0.5
            assignments.append(a)

            q = np.random.rand(edges[i].shape[0], edges[j].shape[0]) - 0.5
            quadratic_terms.append(q)

    instance = mgm.multigraph_matching_input(assignments, quadratic_terms, edges)

    solver = mgm.multigraph_matching_message_passing_solver(['','--multigraphMatchingRoundingMethod','MCF_PS','-v','1','--tighten', '--tightenConstraintsPercentage','0.1','--tightenInterval 50', '--tightenIteration 100', '--tightenReparametrization','uniform:0.5'])
    solver.construct(instance)
    solver.solve()
    result = solver.result()
    result_val = instance.evaluate(result)
    [assignment_masks, quadratic_masks] = result.result_masks(assignments, quadratic_terms, edges)

    result_mask_val = 0.0
    for i in range(len(assignment_masks)):
        result_mask_val = result_mask_val + np.multiply(assignment_masks[i], assignments[i]).sum().sum()
        result_mask_val = result_mask_val + np.multiply(quadratic_masks[i], quadratic_terms[i]).sum().sum()
        
    if(abs(result_val - result_mask_val) > 1e-8):
        raise AssertionError("evaluation difference between native and mask evaluation")


random.seed(1791)

for n in range(3,10):
    nr_nodes = []
    for j in range(n):
        nr_nodes.append(random.randint(2,10))
    construct_mgm_problem(nr_nodes)

    
construct_mgm_problem([3,4,5])

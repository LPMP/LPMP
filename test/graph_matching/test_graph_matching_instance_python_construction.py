import graph_matching_py as gm
import numpy as np


def build_graph(edges_left, edges_right, costs_unary, costs_binary):
    """
        edges_left: list of tuples, orientation matters, no loops [(0,1), (2,1), ...]
        edges_right: the same as above

        costs_unary: np.array(n_1, n_2) where n_i = |V_i|
        cost_binary: np.array(m_1, m_2) where m_i = |E_i|

        return: instance of graph_matching_input
    """
    instance = gm.graph_matching_input()

    # tmp array to store the index of unary assignments
    indices = np.zeros_like(costs_unary, dtype=int)

    it = np.nditer(costs_unary, flags=['multi_index'])
    for cost, n in zip(it, range(costs_unary.size)):
        instance.add_assignment(*it.multi_index, cost)
        indices[it.multi_index] = n
        # print(n, it.multi_index, cost)

    for lei, le in enumerate(edges_left):
        for rei, re in enumerate(edges_right):
            ind1, ind2 = indices[le[0], re[0]], indices[le[1], re[1]]
            instance.add_quadratic_term(ind1, ind2, costs_binary[lei, rei])
            # print(ind1, ind2, costs_2nd[lei, rei])

    return instance


assignments = np.array([[0,1,2,3],
                        [4,5,6,7],
                        [8,9,10,11],
                        [12,13,14,15]], dtype=np.double)

edges_left = np.array([[0,1], [1,2], [2,3]], dtype=np.intc)
edges_right = np.array([[0,1], [1,2], [2,3]], dtype=np.intc)

quadratic_terms = np.array([[-10,-11,-12],
                            [-13,-14,-15],
                            [-16,-17,-18]], dtype=np.double)

instance = gm.graph_matching_input(assignments, quadratic_terms, edges_left, edges_right)
instance_check = build_graph(edges_left, edges_right, assignments, quadratic_terms)

labeling = gm.graph_matching_labeling([0,1,2,3])
[assignment_mask, quadratic_terms_mask] = labeling.result_masks(assignments, quadratic_terms, edges_left, edges_right)

if(instance.evaluate(labeling) != -12):
    raise AssertionError("instance evaluation gives wrong value")
    
if(instance_check.evaluate(labeling) != -12):
    raise AssertionError("instance from batch construction gives wrong value")
if(sum(sum(np.multiply(assignment_mask, assignments))) + sum(sum(np.multiply(quadratic_terms_mask, quadratic_terms))) != -12):
    raise AssertionError("mask solution not correct")






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

def test_graph_construction(edges_left, edges_right, costs_unary, costs_binary):

    nr_left_nodes = costs_unary.shape[0]
    nr_right_nodes = costs_unary.shape[1]
    instance = gm.graph_matching_input(costs_unary, costs_binary, edges_left, edges_right)
    instance_check = build_graph(edges_left, edges_right, assignments, quadratic_terms)

    labeling = gm.graph_matching_labeling(list(range(nr_left_nodes)))
    [assignment_mask, quadratic_terms_mask] = labeling.result_masks(costs_unary, costs_binary, edges_left, edges_right)

    labeling_cost = instance.evaluate(labeling)
    if(instance_check.evaluate(labeling) != labeling_cost):
        raise AssertionError("instance from batch construction gives wrong value")
    if(sum(sum(np.multiply(assignment_mask, assignments))) + sum(sum(np.multiply(quadratic_terms_mask, quadratic_terms))) != labeling_cost):
        raise AssertionError("mask solution not correct")



assignments = np.array([[0,1,2,3],
                        [4,5,6,7],
                        [8,9,10,11],
                        [12,13,14,15]], dtype=np.double)

edges_left = np.array([[0,1], [1,2], [2,3]], dtype=np.intc)
edges_right = np.array([[0,1], [1,2], [2,3]], dtype=np.intc)

quadratic_terms = np.array([[-10,-11,-12],
                            [-13,-14,-15],
                            [-16,-17,-18]], dtype=np.double)

#test_graph_construction(edges_left, edges_right, assignments, quadratic_terms)



edges_left = np.array([[0, 1], [1, 0]], dtype=np.intc)
edges_right = np.array([[0, 2]], dtype=np.intc)
assignments = np.array([[ -5.0, -13.0, -17.0], [ -9.0, -13.0, -15.0]], dtype=np.double)
quadratic_terms = np.array([[-31.0], [-17.0]], dtype=np.double)
test_graph_construction(edges_left, edges_right, assignments, quadratic_terms)

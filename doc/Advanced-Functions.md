# Advanced LPMP capabilities

* [Advanced factor functions](#advanced-factor-functions)
* [Advanced message functions](#advanced-message-functions)
* [Advanced problem constructor functions](#advanced-problem-constructor-functions)
* [Advanced LP functions](#advanced-LP-functions)
* [Lagrange decomposition into trees](#lagrange-decomposition-into-trees)

## Advanced factor functions

In the basic [QPBO factors](Getting-Started#factors) we only provided minimally functioning factors. Additional functions can be provided by the factor and utilized during optimization.

* `void MaximizePotentialAndComputePrimal()`: In some algorithms, whenever a certain factor is visited, a primal solution in that factor is computed by the given function.
* `template<typename ARRAY> void apply(ARRAY& a) const`: When using the proximal bundle method, a subgradient must be computed. The function apply writes 1s into the array `a` in the positions corresponding to the current primal.
* `template<typename EXTERNAL_SOLVER> void construct_constraints(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::${datatype}...) const`: LPMP can export the problem to other ILP solvers and into a text file format. To use this functionality, each factor must provide a function that constructs constraints that correspond to the subproblem the factor solves. The first argument is the external solver, while the latter correspond to the variables exported via the `export_variables` member function.
* `template<typename EXTERNAL_SOLVER> void convert_primal(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::${datatype}...)`: Given a solution from an external solver, extract the primal solution and set it in the factor.

## Advanced message functions

In the basic [QPBO message](Getting-Started#message) we only provided a minimally working message. Additional functions can be provided by the message and utilized during optimization.

* `template<typename LEFT_FACTOR, typename RIGHT_FACTOR> bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)`; Given a primal solution in the left factor, the function sets as many variables in the right factor as possible so that the coupling constraints enforced by the message are fulfilled. Returns true if the primal solution was changed.
* `template<typename LEFT_FACTOR, typename RIGHT_FACTOR> bool ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)`: Same as above but in the other direction.
* `template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR> void construct_constraints(SOLVER& s, const LEFT_FACTOR& l, typename SOLVER::{datatype}... left_variables, const RIGHT_FACTOR& r, typename SOLVER::${datatype}... right_variables) const`: Similarly as for the factor, this function is used to export the constraints that the message represents.
* `template<typename LEFT_FACTOR, typename RIGHT_FACTOR> bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const`: return whether primal solution in left and right factor satisfy the coupling constraints of the respective message.

## Advanced problem constructor functions

In the basic [QPBO problem constructor](Getting-Started#problem-constructor) we only provided a function to construct the problem given the input. More functions are possible which extend the optimization capabilities.

* `void Begin()`: Before optimization starts and after all factors and messages have been added, the `Begin` function is called.
* `void End()`: After optimization terminated, the `End` function is called.
* `std::size_t Tighten()`: Periodically, as indicated by the visitor, the LP-relaxation LPMP optimizes can be tightened and `Tighten` is called to that end.
* `void ComputePrimal()`: Periodically, as governed by the visitor, a primal solution can be decoded by calling the `ComputePrimal()` function.

## Advanced LP functions

We have seen in the basic [QPBO problem constructor](Getting-Started#problem-constructor) how factors and messages can be added to the LP class. LPMP iterates in a forward and backward direction over factors during optimization and the LP class supports indicating in which order this shall happen via the following functions.

* `void add_forward_pass_factor_relation(FactorTypeAdapter* f1, FactorTypeAdapter* f2)`: Indicate that factor `f1` must be processed before factor `f2` in the forward pass.
* `void add_backward_pass_factor_relation(FactorTypeAdapter* f1, FactorTypeAdapter* f2)`: Indicate that factor `f1` must be processed before factor `f2` in the backward pass.
* `void add_factor_relation(FactorTypeAdapter* f1, FactorTypeAdapter* f2)`: Equivalent to `add_forward_pass_factor_relation(f1,f2)` and `add_backward_pass_factor_relation(f2,f1)`.  

## Lagrange decomposition into trees

The standard LPMP decomposition introduces Lagrange multipliers for every message.
 This can lead to a large number of Lagrange multipliers.
 Subgradient based methods are not efficient in such situations, preferring Lagrange decompositions into larger subproblems.
 To this end LPMP offers to build tree-structured subproblems with the `factor_tree` class found [here](/include/tree_decomposition.hxx). 
The `template<typename MESSAGE> void add_message(MESSAGE* msg, Chirality c)` member function add a message to a tree. 
The parameter `Chirality c` can be either set to `Chirality::left` or `Chirality::right` and denotes whether the left or right factor are nearer the root of the tree.
Finally, a collection of trees that cover all messages can be added to and solved by e.g.\ the proximal bundle method [`LP_tree_FWMAP`](/include/LP_FWMAP).


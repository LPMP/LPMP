#ifndef LPMP_LOCAL_PROBLEM_HXX
#define LPMP_LOCAL_PROBLEM_HXX

#include "vector.hxx"
#include <algorithm>
#include <array>

namespace LPMP {


   // graph matching subproblem on a subset of nodes solved with explicit enumeration or branch and bound
   class graph_matching_local_subproblem {
      public:
         template<typename UNARY_ITERATOR, typename ASSIGNMENT_ITERATOR, typename PAIRWISE_ITERATOR>
            graph_matching_local_subproblem(
                  UNARY_ITERATOR unary_size_begin, UNARY_ITERATOR unary_size_end,
                  ASSIGNMENT_ITERATOR assignment_begin, ASSIGNMENT_ITERATOR assignment_end,
                  PAIRWISE_ITERATOR pairwise_idx_begin, PAIRWISE_ITERATOR pairwise_idx_end)
            : unary(unary_size_begin, unary_size_end),
            assignment(unary_size_begin, unary_size_end) 
            {
               const INDEX no_unaries = std::distance(unary_size_begin, unary_size_end);
               assert(3 <= no_unaries && no_unaries <= 4);
               assert(std::distance(assignment_begin, assignment_end) == no_unaries);
               // make initial costs zero
               for(INDEX i=0; i<unary.size(); ++i) {
                  for(INDEX j=0; j<unary[i].size(); ++j) {
                     unary(i,j) = 0.0;
                  }
               }

               // populate assignments
               for(INDEX i=0; i<assignment.size(); ++i) {
                  auto assignments = *(assignment_begin+i);
                  for(INDEX j=0; j<assignment[i].size(); ++j) {
                     assignment(i,j) = assignments[j];
                  }
               }

               no_pairwise_potentials = std::distance(pairwise_idx_begin, pairwise_idx_end);
               pairwise = (pairwise_entry*) malloc(no_pairwise_potentials * sizeof(pairwise_entry));// pairwise_entry[no_pairwise_potentials]; // to do: use stack allocator
               for(INDEX i=0; i<no_pairwise_potentials; ++i) {
                 auto idx = *(pairwise_idx_begin + i);
                 pairwise[i].first_node = idx[0];
                 pairwise[i].second_node = idx[1];
                 new (&pairwise[i].cost) matrix<REAL>(no_labels(idx[0]), no_labels(idx[1]), 0.0);
               }

               assert(no_nodes() <= 4);
               std::array<INDEX,3> pairwise_forward_size({0,0,0});
               std::array<INDEX,3> pairwise_backward_size({0,0,0});
               for(auto idx=pairwise_idx_begin; idx!=pairwise_idx_end; ++idx) {
                  pairwise_forward_size[(*idx)[0]]++; 
                  pairwise_backward_size[(*idx)[1]-1]++;
               }

               pairwise_forward.resize(pairwise_forward_size.begin(), pairwise_forward_size.begin()+no_unaries-1);
               pairwise_backward.resize(pairwise_backward_size.begin(), pairwise_backward_size.begin()+no_unaries-1);

               // populate pointers to pairwise potentials
               std::fill(pairwise_forward_size.begin(), pairwise_forward_size.end(), 0);
               std::fill(pairwise_backward_size.begin(), pairwise_backward_size.end(), 0);
               for(INDEX idx=0; idx<no_pairwise_potentials; ++idx) {
                 const INDEX i = (*(pairwise_idx_begin+idx))[0];
                 const INDEX j = (*(pairwise_idx_begin+idx))[1];
                 pairwise_forward(i,pairwise_forward_size[i]++) = &pairwise[idx];
                 pairwise_backward(j-1,pairwise_backward_size[j]++) = &pairwise[idx];
               } 
            }

         ~graph_matching_local_subproblem()
         {
           if(pairwise != nullptr) { 
             INDEX no_pairwise_potentials = 0;
             for(INDEX i=0; i<pairwise_forward.size(); ++i) {
               no_pairwise_potentials += pairwise_forward[i].size();
             }

             for(INDEX i=0; i<no_pairwise_potentials; ++i) {
               pairwise[i].cost.matrix<REAL>::~matrix();
             } 

             delete[] pairwise; 
           }
         }

         graph_matching_local_subproblem(const graph_matching_local_subproblem& o) { assert(false); } // we would need to do a deep copy because of the pointers in pariwise_2

         // needs to be implemented
         void init_primal() {};
         template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { assert(false); }
         template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) 
         { 
           for(INDEX i=0; i<unary.size(); ++i) {
             ar( binary_data<REAL>(unary[i].begin(), unary[i].size()) ); 
           }
           for(INDEX i=0; i<no_pairwise_potentials; ++i) {
             ar( pairwise[i].cost );
           }
         } 

         auto export_variables() { return std::tie(); }

         template<typename ITERATOR>
         bool feasible(ITERATOR begin, ITERATOR end) const
         {
            assert(std::distance(begin,end) == no_nodes());
            auto it=begin;
            for(INDEX i=0; i<no_nodes(); ++i, ++begin) {
               for(INDEX i=0; i<no_nodes(); ++i) {
                  for(INDEX j=i+1; j<no_nodes(); ++j) {
                     if(assignment(i,*(begin+i)) == assignment(j,*(begin+j))) {
                        return false;
                     }
                  }
               }
            }
            return true;
         }

         REAL EvaluatePrimal() const { return evaluate(primal.begin(), primal.end()); }

         template<typename ITERATOR>
         REAL evaluate(ITERATOR begin, ITERATOR end) const
         {
            REAL cost = unary(0,*begin);
            assert(std::distance(begin,end) == no_nodes());
            if(!feasible(begin,end)) { return std::numeric_limits<REAL>::infinity(); }
            auto it=begin;
            for(INDEX idx=0; idx<no_pairwise_potentials; ++idx) {
              const INDEX i = pairwise[idx].first_node;
              const INDEX j = pairwise[idx].first_node;
              const auto& pot = pairwise[idx].cost;
              cost += pot(*(begin+i), *(begin+j));
            }
            return cost;
         }

         void MaximizePotentialAndComputePrimal()
         {
            lower_bound_exhaustive_search(primal.begin());
         }
         REAL LowerBound() const
         {
           assert(no_nodes() <= 4);
           std::array<INDEX,4> sol_tmp;
           assert(std::abs(lower_bound_exhaustive_search(sol_tmp.begin()) - lower_bound_branch_and_bound(sol_tmp.begin())) <= eps);
           return lower_bound_exhaustive_search(sol_tmp.begin());
         }

         // get any feasible solution by depth first search.
         // Improvement: explore labels according to the order of their unary cost value
         template<typename ITERATOR>
         void feasible_solution(ITERATOR solution_begin) const
         {
           assert(no_nodes() <= 4);
           std::array<SIGNED_INDEX,4> solution;
           std::fill(solution.begin(), solution.begin() + no_nodes(), -1);
            SIGNED_INDEX cur_node=0;
            while(1) {
               if(cur_node < 0) { break; }
               if(cur_node == no_nodes()) { // compare current solution against best one until now
                  std::copy(solution.begin(), solution.end(), solution_begin);
                  assert(feasible(solution.begin(), solution.end()));
                  return;
               } else { // explore search tree
                  solution[cur_node]++;
                  // backtrack if all labels have been visited
                  if(solution[cur_node] == no_labels(cur_node)) {
                     solution[cur_node] = -1;
                     --cur_node;
                     continue;
                  }
                  // backtrack if current partial solution has a conflict
                  for(INDEX i=0; i<cur_node; ++i) {
                     if(assignment(i,solution[i]) == assignment(cur_node,solution[cur_node])) {
                        solution[cur_node] = -1;
                        --cur_node; // conflict detected, backtrack until we can find a node whose label can be incremented 
                        continue;
                     }
                  } 
                  ++cur_node; 
               } 
            } 
            assert(false); // no feasible solution exists at all.
         } 

         template<typename ITERATOR>
         bool forward_ICM(ITERATOR solution_begin) const
         {
            assert(feasible(solution_begin, solution_begin + no_nodes()));
            bool changed = false;
            // forward pass: greedily improve labels 0,1,...
            for(INDEX i=0; i<no_nodes(); ++i) {
               REAL best_cost = std::numeric_limits<REAL>::infinity();
               for(INDEX x=0; x<no_labels(i); ++i) {
                  REAL cost = unary(i,x);
                  if(i < no_nodes()-1) {
                     for(INDEX j=0; j<pairwise_forward[i].size(); ++j) {
                        const INDEX second_node = pairwise_forward(i,j)->second_node;
                        cost += pairwise_forward(i,j)->cost(x, solution_begin[second_node]); 
                     }
                  }
                  if(i > 0) {
                     for(INDEX j=0; j<pairwise_backward[i-1].size(); ++j) {
                        const INDEX first_node = pairwise_backward(i-1,j)->first_node;
                        cost += pairwise_backward(i-1,j)->cost(solution_begin[first_node], x); 
                     }
                  }
                  bool feasible = true;
                  for(INDEX prev_node=0; prev_node<i; ++prev_node) {
                     if(assignment(prev_node, solution_begin[prev_node]) == assignment(i, x)) {
                        feasible = false;
                        break;
                     }
                  } 
                  for(INDEX next_node=i+1; next_node<no_nodes(); ++next_node) {
                     if(assignment(next_node, solution_begin[next_node]) == assignment(i, x)) {
                        feasible = false;
                        break;
                     }
                  } 
                  if(feasible && cost < best_cost) {
                     best_cost = cost;
                     if(x != solution_begin[i]) { changed = true; }
                     solution_begin[i] = x;
                  }
               } 
            }
            assert(feasible(solution_begin, solution_begin + no_nodes()));
            return changed;
         }
         template<typename ITERATOR>
         bool backward_ICM(ITERATOR solution_begin) const
         {
            assert(feasible(solution_begin, solution_begin + no_nodes()));
            bool changed = false;
            // backward pass: greedily improve labels no_nodes()-1,no_nodes()-2,...
            for(SIGNED_INDEX i=no_nodes()-1; i>=0; --i) {
               REAL best_cost = std::numeric_limits<REAL>::infinity();
               for(INDEX x=0; x<no_labels(i); ++i) {
                  REAL cost = unary(i,x);
                  if(i < no_nodes()-1) {
                    assert(false); // the order is the same as for forward_ICM
                     for(INDEX j=0; j<pairwise_forward[i].size(); ++j) {
                        const INDEX second_node = pairwise_forward(i,j)->second_node;
                        cost += pairwise_forward(i,j)->cost(x, solution_begin[second_node]); 
                     }
                  }
                  if(i > 0) {
                     for(INDEX j=0; j<pairwise_backward[i-1].size(); ++j) {
                        const INDEX first_node = pairwise_backward(i-1,j)->first_node;
                        cost += pairwise_backward(i-1,j)->cost(solution_begin[first_node], x); 
                     }
                  }
                  bool feasible = true;
                  for(INDEX prev_node=0; prev_node<i; ++prev_node) {
                     if(assignment(prev_node, solution_begin[prev_node]) == assignment(i, x)) {
                        feasible = false;
                        break;
                     }
                  } 
                  for(INDEX next_node=i+1; next_node<no_nodes(); ++next_node) {
                     if(assignment(next_node, solution_begin[next_node]) == assignment(i, x)) {
                        feasible = false;
                        break;
                     }
                  }
                  if(feasible && cost < best_cost) {
                     best_cost = cost;
                     if(x != solution_begin[i]) { changed = true; }
                     solution_begin[i] = x;
                  }
               } 
            }
            assert(feasible(solution_begin, solution_begin + no_nodes()));
            return changed;
         }

         using lb_data_type = two_dim_variable_array<REAL>;

         template<typename ITERATOR>
         REAL lb(const lb_data_type& lb, ITERATOR sol_begin, ITERATOR sol_end) const
         {
            REAL val = 0.0;
            const INDEX d = std::distance(sol_begin, sol_end);
            assert(d < no_nodes());
            for(INDEX i=d+1; i<no_nodes(); ++i) {
               REAL best_cost = std::numeric_limits<REAL>::infinity();
               for(INDEX x=0; x<no_labels(i); ++x) { // possibly exclude labels that lead to infeasibility
                  REAL cur_node_cost = 0.0;
                  // unaries and pairwise potentials that have support nodes with index >= d
                  cur_node_cost += lb(i,x);
                  // pairwise potentials where first variable is already chosen
                  for(INDEX j=0; j<pairwise_backward[i-1].size(); ++j) {
                     const INDEX first_node = pairwise_backward(i-1,j)->first_node;
                     if(first_node < d) {
                        cur_node_cost += pairwise_backward(i-1,j)->cost( *(sol_begin+first_node), x);
                     }
                  }
                  best_cost = std::min(best_cost, cur_node_cost);
               }
               val += best_cost;
            }
            return val;
         }
         lb_data_type compute_lb_data() const
         {
            // marginalize out the second variable in all pairwise potentials
            lb_data_type min_marg(unary);
            for(INDEX i=0; i<no_nodes()-1; ++i) {
               for(INDEX j=0; j<pairwise_forward[i].size(); ++j) {
                  for(INDEX x=0; x<no_labels(i); ++x) {
                     min_marg(i,x) += pairwise_forward(i,j)->cost.col_min(x);
                  } 
               } 
            } 
            return min_marg;
         }

         template<typename ITERATOR>
         REAL lower_bound_branch_and_bound(ITERATOR solution_begin) const
         {
            vector<INDEX> best_solution(no_nodes());
            // get a feasible solution
            REAL best_cost = evaluate(primal.begin(), primal.end());
            if(best_cost == std::numeric_limits<REAL>::infinity()) {
               feasible_solution(best_solution.begin());
               best_cost = evaluate(best_solution.begin(), best_solution.end());
            } else {
               best_solution = primal;
            }

            // improve feasible solution via ICM
            while(forward_ICM(best_solution.begin()) && backward_ICM(best_solution.begin())) {} 

            // branch and bound
            vector<SIGNED_INDEX> current_solution(no_nodes(), -1);
            vector<REAL> current_cost(no_nodes()+1);
            auto lb_data = compute_lb_data();

            SIGNED_INDEX cur_node=0;
            current_cost[0] = 0.0;
            while(1) {
               if(cur_node < 0) { break; }
               if(cur_node == no_nodes()) { // compare current solution against best one until now
                  if(current_cost.back() < best_cost) {
                     std::copy(current_solution.begin(), current_solution.end(), best_solution.begin());
                     best_cost = current_cost.back();
                  }
                  --cur_node; 
               } else { // explore search tree
                  // backtrack if current partial solution has a conflict
                  for(INDEX i=0; i<cur_node; ++i) {
                     if(assignment(i,current_solution[i]) == assignment(cur_node,current_solution[cur_node])) {
                        current_solution[cur_node] = -1;
                        --cur_node; // conflict detected, backtrack until we can find a node whose label can be incremented 
                        continue;
                     }
                  } 
                  // estimate lower bound
                  const REAL lb_val = current_cost[cur_node+1] + lb(lb_data, current_solution.begin(), current_solution.begin()+cur_node+1);
                  if(best_cost <= lb_val) { // backtrack, no need to explore this part of search tree further
                     current_solution[cur_node]++;
                     if(current_solution[cur_node] == no_labels(cur_node)) {
                        current_solution[cur_node] = -1;
                        --cur_node;
                     } 
                  } else {
                     ++cur_node;
                  }
               }
            }

            std::copy(best_solution.begin(), best_solution.end(), solution_begin);
            return best_cost; 
         }

         template<typename ITERATOR>
         REAL lower_bound_exhaustive_search(ITERATOR solution_begin) const
         {
            vector<SIGNED_INDEX> best_solution(no_nodes());
            REAL best_cost = std::numeric_limits<REAL>::infinity();
            vector<SIGNED_INDEX> current_solution(no_nodes(), -1);
            vector<REAL> current_cost(no_nodes()+1);

            SIGNED_INDEX cur_node=0;
            current_cost[0] = 0.0;
            while(1) {
               if(cur_node < 0) { break; }
               if(cur_node == no_nodes()) { // compare current solution against best one until now
                  if(current_cost.back() < best_cost) {
                     best_solution = current_solution;
                     best_cost = current_cost.back();
                  }
                  --cur_node; 
               } else { // explore search tree
                  current_solution[cur_node]++;
                  // backtrack if all labels have been visited
                  if(current_solution[cur_node] == no_labels(cur_node)) {
                     current_solution[cur_node] = -1;
                     --cur_node;
                     continue;
                  }
                  // backtrack if current partial solution has a conflict
                  for(INDEX i=0; i<cur_node; ++i) {
                     if(assignment(i,current_solution[i]) == assignment(cur_node,current_solution[cur_node])) {
                        current_solution[cur_node] = -1;
                        --cur_node; // conflict detected, backtrack until we can find a node whose label can be incremented 
                        continue;
                     }
                  } 
                  // record cost
                  current_cost[cur_node+1] = current_cost[cur_node];
                  for(INDEX j=0; j<pairwise_backward[cur_node-1].size(); ++j) {
                     const INDEX first_node = pairwise_backward(cur_node-1,j)->first_node;
                     current_cost[cur_node+1] += pairwise_backward(cur_node-1,j)->cost(current_solution[first_node], current_solution[cur_node]); 
                  }
                  ++cur_node; 
               } 
            } 

            std::copy(best_solution.begin(), best_solution.end(), solution_begin);
            return best_cost;
         } 

         INDEX no_nodes() const { return unary.size(); }
         INDEX no_labels(const INDEX i) const { assert(i < unary.size()); return unary[i].size(); }

         template<typename SOLVER>
         void construct_constraints(SOLVER& s)
         {}

         template<typename SOLVER>
         void convert_primal(SOLVER& s)
         {}

         vector<INDEX> primal;

         two_dim_variable_array<REAL> unary;
         // for each node except the first one, a list of pairwise potentials is stored with first_node denoting its first.
         struct pairwise_entry { INDEX first_node; INDEX second_node; matrix<REAL> cost; };
         pairwise_entry* pairwise = nullptr;
         //two_dim_variable_array<pairwise_entry> pairwise;
      private:

         // the same but reversed
         //struct pairwise_entry_2 { INDEX second_node; matrix<REAL>* cost; };
         //two_dim_variable_array<pairwise_entry_2> pairwise_forward;

         two_dim_variable_array<INDEX> assignment;
         // possibly also store a vector with compressed assignment information such that taken assignments can be looked up fast

         //struct pairwise_entry {
         //  std::array<INDEX,2> idx;
         //  matrix<REAL> cost;
         //}

         //pairwise_entry* pairwise;
         two_dim_variable_array<pairwise_entry*> pairwise_forward; 
         two_dim_variable_array<pairwise_entry*> pairwise_backward; 
         INDEX no_pairwise_potentials;
   };

   // left is unary factor, right is local problem
   class unary_local_problem_message {
      public:
         unary_local_problem_message(const INDEX _unary_no) : unary_no(_unary_no) {}

         template<typename LEFT_POT, typename MSG_ARRAY>
            void send_message_to_right(const LEFT_POT& leftPot, MSG_ARRAY& msg, const REAL omega) 
            {
               for(INDEX i=0; i<leftPot.size(); ++i) {
                  msg[i] -= omega*(leftPot[i]);
               }
            }

         // currently this function must be present, although SendMessagesToLeft is called
         template<typename RIGHT_POT, typename MSG_ARRAY>
            void send_message_to_left(const RIGHT_POT& r, MSG_ARRAY& msg, const REAL omega) 
            {
               assert(false);
            }

         template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
            void ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
            {
               l.primal() = r.primal[unary_no];
            }

         template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
            void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
            {
               assert(false);
            }

         template<typename LEFT_FACTOR>
            void RepamLeft(LEFT_FACTOR& left, const REAL msg, const INDEX dim) 
            {
               assert(dim < left.size());
               left[dim] += msg;
            }

         template<typename RIGHT_FACTOR>
            void RepamRight(RIGHT_FACTOR& right, const REAL msg, const INDEX dim) 
            {
               right.unary(unary_no,dim) += msg;
            }

         template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
         void construct_constraints(SOLVER& s, LEFT_FACTOR& unary, typename SOLVER::vector unary_variables, RIGHT_FACTOR& r)
         {}

      private:
         const INDEX unary_no;
   };

   // left is pairwise factor, right is local problem
   class pairwise_local_problem_message {
      public:
         pairwise_local_problem_message(const INDEX unary1, const INDEX unary2, const graph_matching_local_subproblem& f) 
         {
            pairwise_idx = std::numeric_limits<INDEX>::max();
            for(INDEX i=0; ; ++i) {
              if(f.pairwise[i].first_node == unary1 && f.pairwise[i].second_node == unary2) {
                pairwise_idx = i;
                break;
              }
            }
            assert(pairwise_idx < std::numeric_limits<INDEX>::max()); 
         }

         template<typename LEFT_POT, typename MSG_ARRAY>
            void send_message_to_right(const LEFT_POT& l, MSG_ARRAY& msg, const REAL omega) 
            {
               msg -= omega*l;
            }

         template<typename RIGHT_POT, typename MSG_ARRAY>
            void send_message_to_left(const RIGHT_POT& r, MSG_ARRAY& msg, const REAL omega) 
            {
               assert(false);
            }

         template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
            void ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
            {
               const INDEX unary1 = r.pairwise[pairwise_idx].first_node;
               const INDEX unary2 = r.pairwise[pairwise_idx].second_node;
               l.primal()[0] = r.primal[unary1];
               l.primal()[1] = r.primal[unary2];
            }

         template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
            void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
            {
               assert(false);
            }


         template<typename LEFT_FACTOR, typename MSG>
            void RepamLeft(LEFT_FACTOR& l, const MSG& msgs)
            {
               // do zrobienia: possibly use counter
               for(INDEX x1=0; x1<l.dim1(); ++x1) {
                  for(INDEX x2=0; x2<l.dim2(); ++x2) {
                     l.cost(x1,x2) += normalize( msgs(x1,x2) );
                     assert(!std::isnan(l(x1,x2)));
                  }
               }
            }

         template<typename RIGHT_FACTOR, typename MSG>
            void RepamRight(RIGHT_FACTOR& r, const MSG& msgs)
            {
               auto& cost = r.pairwise[pairwise_idx].cost;
               assert(msgs.dim1() == cost.dim1() && msgs.dim2() == cost.dim2());
               for(INDEX x1=0; x1<cost.dim1(); ++x1) {
                  for(INDEX x2=0; x2<cost.dim2(); ++x2) {
                     cost(x1,x2) += normalize( msgs(x1,x2) );
                     assert(!std::isnan(cost(x1,x2)));
                  }
               }
            }

         template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
         void construct_constraints(SOLVER& s, 
             LEFT_FACTOR& pairwise, typename SOLVER::vector left_unary_variables, typename SOLVER::vector right_unary_variables, typename SOLVER::vector pairwise_variables,
             RIGHT_FACTOR& r)
         {}

      private:
         INDEX pairwise_idx;
   };

   template<typename LOCAL_SUBPROBLEM_CONTAINER, typename UNARY_LOCAL_SUBPROBLEM_MESSAGE, typename PAIRWISE_LOCAL_SUBPROBLEM_MESSAGE>
   class local_subproblem_constructor {
      public:
      using FMC = typename LOCAL_SUBPROBLEM_CONTAINER::FMC;

         template<typename SOLVER>
            local_subproblem_constructor(SOLVER& solver) : lp_(&solver.GetLP()) {}

         template<typename ITERATOR>
         bool has_factor(ITERATOR begin, ITERATOR end)
         {
            const INDEX n = std::distance(begin, end);
            assert(3 <= n && n <= 4);
            assert(std::is_sorted(begin,end));
            if(n == 3) {
               std::array<INDEX,3> idx({*begin,*(begin+1),*(begin+2)});
               return local_subproblems3.find(idx) != local_subproblems3.end();
            } else {
               std::array<INDEX,4> idx({*begin,*(begin+1),*(begin+2),*(begin+3)});
               return local_subproblems4.find(idx) != local_subproblems4.end();
            } 
         }
         template<typename ITERATOR, typename MRF_CONSTRUCTOR>
         LOCAL_SUBPROBLEM_CONTAINER* add_local_subproblem(ITERATOR begin, ITERATOR end, const MRF_CONSTRUCTOR& mrf, factor_tree<FMC>* t = nullptr)
         {
            assert(!has_factor(begin,end));

            const INDEX n = std::distance(begin, end);
            assert(3 <= n && n <= 4);
            assert(std::is_sorted(begin,end));
            // collect all unary sizes
            std::array<INDEX,4> size;
            for(INDEX i=0; i<n; ++i) {
               size[i] = mrf.GetNumberOfLabels(*(begin+i));
            }

            // collect all pairwise potentials
            std::vector<std::array<INDEX,2>> pairwise_potentials; // first variable, second variable
            for(INDEX i=0; i<n; ++i) {
               for(INDEX j=i+1; j<n; ++j) {
                  if(mrf.HasPairwiseFactor(*(begin+i), *(begin+j))) {
                     pairwise_potentials.push_back({i,j}); 
                  }
               }
            }

            // get assignments variables for each node
            std::array<std::vector<INDEX>,4> assignments;
            for(INDEX i=0; i<n; ++i) {
               assignments[i] = mrf.assignment_nodes(*(begin+i));
            }

            auto* f = lp_->template add_factor<LOCAL_SUBPROBLEM_CONTAINER>(size.begin(), size.begin()+n, assignments.begin(), assignments.begin()+n, pairwise_potentials.begin(), pairwise_potentials.end());

            // connect with unary potentials
            for(INDEX i=0; i<n; ++i) {
               auto* u = mrf.GetUnaryFactor(*(begin+i));
               auto* m = lp_->template add_message<UNARY_LOCAL_SUBPROBLEM_MESSAGE>(u, f, i); 
               if(t != nullptr) { t->add_message(m, Chirality::right); }
            }

            // connect with pairwise potentials
            for(auto idx : pairwise_potentials) {
               auto* p = mrf.GetPairwiseFactor(*(begin+idx[0]), *(begin+idx[1]));
               auto* m = lp_->template add_message<PAIRWISE_LOCAL_SUBPROBLEM_MESSAGE>(p, f, pairwise_local_problem_message(idx[0], idx[1], *f->get_factor()));
               if(t != nullptr) { t->add_message(m, Chirality::right); }
            }

            if(n == 3) {
               std::array<INDEX,3> idx({*begin,*(begin+1),*(begin+2)});
               local_subproblems3.insert(std::make_pair(idx,f));
            } else {
               assert(n == 4);
               std::array<INDEX,4> idx({*begin,*(begin+1),*(begin+2),*(begin+3)});
               local_subproblems4.insert(std::make_pair(idx,f));
               // also add to local_subproblems3 all subsets of cardinality 3, if they are not there already
               {
                  std::array<INDEX,3> sub_idx({*(begin+1),*(begin+2),*(begin+3)});
                  if(!has_factor(sub_idx.begin(), sub_idx.end())) { local_subproblems3.insert(std::make_pair(sub_idx,f)); }
               }
               {
                  std::array<INDEX,3> sub_idx({*begin,*(begin+2),*(begin+3)});
                  if(!has_factor(sub_idx.begin(), sub_idx.end())) { local_subproblems3.insert(std::make_pair(sub_idx,f)); }
               }
               {
                  std::array<INDEX,3> sub_idx({*begin,*(begin+1),*(begin+3)});
                  if(!has_factor(sub_idx.begin(), sub_idx.end())) { local_subproblems3.insert(std::make_pair(sub_idx,f)); }
               }
               {
                  std::array<INDEX,3> sub_idx({*begin,*(begin+1),*(begin+2)});
                  if(!has_factor(sub_idx.begin(), sub_idx.end())) { local_subproblems3.insert(std::make_pair(sub_idx,f)); }
               }
            }

            return f;
         }

         template<typename MRF_CONSTRUCTOR>
         std::vector<factor_tree<FMC>> add_local_subproblems(const MRF_CONSTRUCTOR& mrf)
         {
            // iterate over all 4-cycles and add subproblems

            // iterate over 3-cliques and add subproblems if no subproblem covers it yet.
            return add_local_subproblems_triplet(mrf);
         }

         template<typename MRF_CONSTRUCTOR>
         std::vector<factor_tree<FMC>> add_local_subproblems_triplet(const MRF_CONSTRUCTOR& mrf, factor_tree<FMC>* t = nullptr)
         {
            vector<INDEX> adjacency_list_count(mrf.GetNumberOfVariables(),0);
            // first determine size for adjacency_list
            for(INDEX i=0; i<mrf.GetNumberOfPairwiseFactors(); ++i) {
               auto idx = mrf.GetPairwiseVariables(i);
               adjacency_list_count[idx[0]]++;
               adjacency_list_count[idx[1]]++;
            }
            two_dim_variable_array<INDEX> adjacency_list(adjacency_list_count.begin(), adjacency_list_count.end());
            std::fill(adjacency_list_count.begin(), adjacency_list_count.end(), 0);
            for(INDEX i=0; i<mrf.GetNumberOfPairwiseFactors(); ++i) {
               auto edge_idx = mrf.GetPairwiseVariables(i);
               adjacency_list[edge_idx[0]][adjacency_list_count[edge_idx[0]]] = edge_idx[1];
               adjacency_list_count[edge_idx[0]]++;
               adjacency_list[edge_idx[1]][adjacency_list_count[edge_idx[1]]] = edge_idx[0];
               adjacency_list_count[edge_idx[1]]++;
            }

            // Sort the adjacency list, for fast intersections later 
#pragma omp parallel for schedule(guided)
            for(int i=0; i < adjacency_list.size(); i++) {
               std::sort(adjacency_list[i].begin(), adjacency_list[i].end());
            }

            // Iterate over all of the edge intersection sets
            // do zrobienia: parallelize
            std::vector<factor_tree<FMC>> trees;
            {
               std::vector<INDEX> commonNodes(mrf.GetNumberOfVariables());
               std::vector<triplet_candidate> triplet_candidates_per_thread;
               for(INDEX i=0; i<mrf.GetNumberOfPairwiseFactors(); ++i) {
                  auto edge_idx = mrf.GetPairwiseVariables(i);

                     // Now find all neighbors of both i and j to see where the triangles are
                     // TEMP TEMP -- fails at i=0, j=1, on i==3.
                     auto intersects_iter_end = std::set_intersection(
                           adjacency_list[edge_idx[0]].begin(), adjacency_list[edge_idx[0]].end(),
                           adjacency_list[edge_idx[1]].begin(), adjacency_list[edge_idx[1]].end(),
                           commonNodes.begin());

                  for(auto n=commonNodes.begin(); n != intersects_iter_end; ++n) {
                     const INDEX k = *n;

                     // Since a triplet shows up three times as an edge plus
                     // a node, we only consider it for the case when i<j<k 
                     if(!(edge_idx[1]<k)) { continue; }

                     std::array<INDEX,3> triplet_idx({edge_idx[0], edge_idx[1], k});
                     if(!has_factor(triplet_idx.begin(), triplet_idx.end())) {
                       factor_tree<FMC> t;
                       add_local_subproblem(triplet_idx.begin(), triplet_idx.end(), mrf, &t);
                       trees.push_back(t);
                     }
                  } 
               }
            }
            return trees;
         }


      private:
         std::unordered_map<std::array<INDEX,3>, LOCAL_SUBPROBLEM_CONTAINER*> local_subproblems3;
         std::unordered_map<std::array<INDEX,4>, LOCAL_SUBPROBLEM_CONTAINER*> local_subproblems4; 

         LP<FMC>* lp_;
   };


} // namespace LPMP

#endif // LPMP_LOCAL_PROBLEM_HXX

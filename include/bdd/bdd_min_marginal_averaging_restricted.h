#pragma once

#include "bdd_min_marginal_averaging.h"

namespace LPMP {

    // do message passing on subsets of relevant variables only
    class bdd_min_marginal_averaging_restricted : public bdd_min_marginal_averaging {
        public:
            template<typename ITERATOR>
                double variable_score(ITERATOR marginals_begin, ITERATOR marginals_end, const double th) const;

            template<typename VAR_ITERATOR, typename BDD_MASK>
                void min_marginal_averaging_forward_restricted(VAR_ITERATOR var_begin, VAR_ITERATOR var_end, BDD_MASK bdd_mask);

            template<typename VAR_ITERATOR, typename BDD_MASK>
                void min_marginal_averaging_backward_restricted(VAR_ITERATOR var_begin, VAR_ITERATOR var_end, BDD_MASK bdd_mask);

            void min_marginal_averaging_iteration_restricted();

            template<typename VARIABLES, typename BDD_MASK>
                void min_marginal_averaging_iteration_restricted(VARIABLES variables, BDD_MASK bdd_mask);

        private:
        template<typename ITERATOR>
            std::tuple<two_dim_variable_array<char>,std::vector<std::size_t>> compute_bdd_mask(ITERATOR variable_begin, ITERATOR variable_end) const;
    }; 


   template<typename ITERATOR>
       double bdd_min_marginal_averaging_restricted::variable_score(ITERATOR marginals_begin, ITERATOR marginals_end, const double th) const
        {
            double largest_positive = -std::numeric_limits<double>::infinity();
            double smallest_negative = std::numeric_limits<double>::infinity();
            std::size_t zero_min_marginal_diff;
            for(auto marginals_it = marginals_begin; marginals_it != marginals_end; ++marginals_it) {
                const double marginal_diff = (*marginals_it)[1] - (*marginals_it)[0];
                if(marginal_diff > th)
                    largest_positive = std::max(largest_positive, marginal_diff);
                else if(marginal_diff < -th)
                    smallest_negative = std::min(smallest_negative,marginal_diff); 
                else 
                    zero_min_marginal_diff++;
            }

            if(std::isfinite(largest_positive) && std::isfinite(smallest_negative)) {
                const double lower_bound_gain = std::min(largest_positive, -smallest_negative);
                return lower_bound_gain;
            } else if(zero_min_marginal_diff == std::distance(marginals_begin, marginals_end)) {
                return std::numeric_limits<double>::infinity();
            }
            return 0.0;
        }

    template<typename VAR_ITERATOR, typename BDD_MASK>
        void bdd_min_marginal_averaging_restricted::min_marginal_averaging_forward_restricted(VAR_ITERATOR var_begin, VAR_ITERATOR var_end, BDD_MASK bdd_mask)
        {
            assert(std::is_sorted(var_begin, var_end));
            std::vector<std::array<double,2>> min_marginals;
            for(auto var_it=var_begin; var_it!=var_end; ++var_it) {
                const std::size_t var = *var_it;
                //std::cout << "variable = " << var << "; ";

                // collect min marginals
                min_marginals.clear();
                for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                    if(bdd_mask(var,bdd_index)) {
                        forward_step(var,bdd_index);
                        //std::cout << bdd_index << ", ";
                        min_marginals.push_back(min_marginal(var,bdd_index));
                    }
                }
                //std::cout << "\n";

                // set marginals in each bdd so min marginals match each other
                if(min_marginals.size() > 1) {
                    const std::array<double,2> average_marginal = average_marginals(min_marginals.begin(), min_marginals.end());
                    std::size_t min_marginal_counter = 0;
                    for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                        if(bdd_mask(var,bdd_index)) {
                            set_marginal(var,bdd_index,average_marginal,min_marginals[min_marginal_counter++]);
                        }
                    }
                }
            }
        }

    template<typename VAR_ITERATOR, typename BDD_MASK>
        void bdd_min_marginal_averaging_restricted::min_marginal_averaging_backward_restricted(VAR_ITERATOR var_begin, VAR_ITERATOR var_end, BDD_MASK bdd_mask)
        {
            assert(std::is_sorted(var_begin, var_end, std::greater_equal<std::size_t>()));
            std::vector<std::array<double,2>> min_marginals;

            for(auto var_it=var_begin; var_it!=var_end; ++var_it) {
                const std::size_t var = *var_it;

                // collect min marginals
                min_marginals.clear();
                for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                    if(bdd_mask(var, bdd_index))
                        min_marginals.push_back(min_marginal(var,bdd_index));
                }

                const std::array<double,2> average_marginal = min_marginals.size() > 0 ? average_marginals(min_marginals.begin(), min_marginals.end()) : std::array<double,2>{0.0,0.0};

                // set marginals in each bdd so min marginals match each other
                std::size_t min_marginal_counter = 0;
                for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                    if(bdd_mask(var, bdd_index)) {
                        set_marginal(var,bdd_index,average_marginal,min_marginals[min_marginal_counter++]);
                        backward_step(var, bdd_index);
                    }
                }
            }
        }

    void bdd_min_marginal_averaging_restricted::min_marginal_averaging_iteration_restricted()
    {
        constexpr double th = 1e-4;
        std::vector<std::size_t> variables;

        // do min-marginal averaging iteration and record score of min-marginals
        std::vector<std::array<double,2>> min_marginals;
        std::cout << "nr vars = " << nr_variables() << "\n";
        for(std::size_t var=0; var<nr_variables(); ++var) {

            // collect min marginals
            min_marginals.clear();
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                forward_step(var,bdd_index);
                min_marginals.push_back(this->min_marginal(var,bdd_index)); 
            }
            const double s = variable_score(min_marginals.begin(), min_marginals.end(), th);
            if(s >= th)
                variables.push_back(var);

            const std::array<double,2> average_marginal = average_marginals(min_marginals.begin(), min_marginals.end());

            // set marginals in each bdd so min marginals match each other
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                set_marginal(var,bdd_index,average_marginal,min_marginals[bdd_index]);
            } 
        } 

        std::cout << "#variables with high scores for restricted iterations = " << variables.size() << "\n";
        std::cout << "total # variables = " << nr_variables() << "\n";
        min_marginal_averaging_backward();

        lower_bound_ = compute_lower_bound();
        std::cout << "lower bound after outer iteration = " << lower_bound_ << "\n";

        const auto [bdd_mask, affected_variables] = compute_bdd_mask(variables.begin(), variables.end());
        std::cout << "#variables in restricted iteration = " << affected_variables.size() << "\n";
        for(std::size_t iter=0; iter<10; ++iter) {
            min_marginal_averaging_iteration_restricted(affected_variables, bdd_mask);
        }
        lower_bound_ = compute_lower_bound();
        std::cout << "lower bound after inner iterations = " << lower_bound_ << "\n";
    }

    template<typename VARIABLES, typename BDD_MASK>
        void bdd_min_marginal_averaging_restricted::min_marginal_averaging_iteration_restricted(VARIABLES variables, BDD_MASK bdd_mask)
        {
            min_marginal_averaging_forward_restricted(variables.begin(), variables.end(), bdd_mask);
            min_marginal_averaging_backward_restricted(variables.rbegin(), variables.rend(), bdd_mask); 
        }


    // given variables, mark all associated bdds
    template<typename ITERATOR>
        std::tuple<two_dim_variable_array<char>,std::vector<std::size_t>> bdd_min_marginal_averaging_restricted::compute_bdd_mask(ITERATOR variable_begin, ITERATOR variable_end) const
        {
            two_dim_variable_array<char> marked_bdds(bdd_variables_);
            std::fill(marked_bdds.data().begin(), marked_bdds.data().begin(), 0);
            std::queue<std::size_t> q;
            for(auto variable_it=variable_begin; variable_it!=variable_end; ++variable_it) {
                const std::size_t var = *variable_it;
                for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                    q.push(marked_bdds.index(var,bdd_index));
                    marked_bdds(var,bdd_index) = 1;
                } 
            }

            while(!q.empty()) {
                const std::size_t idx = q.front();
                q.pop();
                assert(marked_bdds.data()[idx] == 1);
                const auto* prev_bdd_variable = bdd_variables_.data()[idx].prev;
                if(prev_bdd_variable != nullptr) {
                    const std::size_t prev_index = std::distance(&bdd_variables_.data()[0], prev_bdd_variable);
                    if(marked_bdds.data()[prev_index] == 0) {
                        q.push(prev_index);
                        marked_bdds.data()[prev_index] = 1;
                    }
                }

                const auto* next_bdd_variable = bdd_variables_.data()[idx].next;
                if(next_bdd_variable != nullptr) {
                    const std::size_t next_index = std::distance(&bdd_variables_.data()[0], next_bdd_variable);
                    if(marked_bdds.data()[next_index] == 0) {
                        q.push(next_index);
                        marked_bdds.data()[next_index] = 1;
                    }
                }
            }

            std::vector<std::size_t> vars;
            for(std::size_t v=0; v<nr_variables(); ++v) {
                for(std::size_t bdd_index=0; bdd_index<nr_bdds(v); ++bdd_index) {
                    if(marked_bdds(v,bdd_index)) {
                        vars.push_back(v);
                        break;
                    }
                }
            }

            assert(std::is_sorted(vars.begin(), vars.end()));
            for(std::size_t i=1; i<vars.size(); ++i) {
                assert(vars[i-1] < vars[i]);
            }

            return {marked_bdds, vars};
        }
}

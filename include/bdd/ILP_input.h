#pragma once

#include "config.hxx"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <cassert>
#include "two_dimensional_variable_array.hxx"
#include "cuthill-mckee.h"
#include "bfs_ordering.hxx"

namespace LPMP {

    class ILP_input {
        public:
        struct weighted_variable {
            int coefficient;
            std::size_t var;
            bool operator<(const weighted_variable& o) const { return var < o.var; }
        };

        struct linear_constraint {
            std::vector<weighted_variable> variables;
            inequality_type ineq;
            int right_hand_side;
            void normalize() { std::sort(variables.begin(), variables.end()); }
        };

        bool var_exists(const std::string& var) const
        {
            return var_name_to_index_.count(var) > 0;
        }

        std::size_t get_var_index(const std::string& var) const
        {
            assert(var_exists(var));
            return var_name_to_index_.find(var)->second;
        }

        std::size_t add_new_variable(const std::string& var)
        {
            assert(!var_exists(var));
            const std::size_t var_index = var_name_to_index_.size();
            var_name_to_index_.insert({var, var_index});
            assert(var_index_to_name_.size() == var_index);
            var_index_to_name_.push_back(var);
            return var_index;
        }

        std::size_t get_or_create_variable_index(const std::string& var)
        {
            if(var_exists(var))
                return get_var_index(var);
            else 
                return add_new_variable(var); 
        }

        std::size_t nr_variables() const
        {
            return var_name_to_index_.size();
        }

        void add_to_objective(const double coefficient, const std::string& var)
        {
            add_to_objective(coefficient, get_or_create_variable_index(var));
        }

        void add_to_objective(const double coefficient, const std::size_t var)
        {
            if(objective_.size() <= var)
                objective_.resize(var+1,0.0);
            objective_[var] += coefficient;
        }

        const std::vector<double>& objective() const { return objective_; }

        double objective(const std::size_t var) const
        {
            if(var >= nr_variables())
                throw std::runtime_error("variable not present");
            if(var >= objective_.size())
                return 0.0;
            return objective_[var];
        }

        double objective(const std::string& var) const
        {
            return objective(get_var_index(var));
        }

        void begin_new_inequality()
        {
            linear_constraints_.push_back({});
        }

        void set_inequality_type(const inequality_type ineq)
        {
            assert(linear_constraints_.size() > 0);
            linear_constraints_.back().ineq = ineq;
        }

        void add_to_constraint(const int coefficient, const std::size_t var)
        {
            assert(linear_constraints_.size() > 0);
            auto& constr = linear_constraints_.back();
            constr.variables.push_back({coefficient, var}); 

            if(constr.variables.size() > 1)
                if(constr.variables.back() < constr.variables[constr.variables.size()-2])
                    constr.normalize();
        }
        void add_to_constraint(const int coefficient, const std::string& var)
        {
            add_to_constraint(coefficient, get_or_create_variable_index(var));
        }

        void set_right_hand_side(const int x)
        {
            assert(linear_constraints_.size() > 0);
            linear_constraints_.back().right_hand_side = x;
        } 

        std::size_t nr_constraints() const
        {
            return linear_constraints_.size();
        }
            
        const auto& constraints() const
        {
            return linear_constraints_;
        }

        template<typename STREAM>
            void write(STREAM& s) const
            {
                s << "Minimize\n";
                for(const auto o : var_name_to_index_) {
                    s << (objective(o.second) < 0.0 ? "- " : "+ ") <<  std::abs(objective(o.second)) << " " << o.first << "\n"; 
                }
                s << "Subject To\n";
                for(const auto& ineq : constraints()) {
                    for(const auto term : ineq.variables) {
                        s << (term.coefficient < 0.0 ? "- " : "+ ") <<  std::abs(term.coefficient) << " " << var_index_to_name_[term.var] << " "; 
                    }

                    switch(ineq.ineq) {
                        case inequality_type::smaller_equal:
                            s << " <= ";
                            break;
                        case inequality_type::greater_equal:
                            s << " >= ";
                            break;
                        case inequality_type::equal:
                            s << " = ";
                            break;
                        default:
                            throw std::runtime_error("inequality type not supported");
                            break;
                    }
                    s << ineq.right_hand_side << "\n";
                }
                s << "Bounds\n";
                s << "Binaries\n";
                for(const auto& v : var_index_to_name_)
                    s << v << "\n";
                s << "End\n";

            }

        private:
            std::vector<linear_constraint> linear_constraints_;
            std::vector<double> objective_;
            std::vector<std::string> var_index_to_name_;
            std::unordered_map<std::string, std::size_t> var_name_to_index_;
        public:
            void reorder_bfs();
            void reorder_Cuthill_McKee(); 

        private:
            two_dim_variable_array<std::size_t> variable_adjacency_matrix() const;
            void reorder(const std::vector<std::size_t>& new_order);
    };

    inline two_dim_variable_array<std::size_t> ILP_input::variable_adjacency_matrix() const
    {
        std::unordered_set<std::array<std::size_t,2>> adjacent_vars;
        for(const auto& l : this->linear_constraints_) {
            for(std::size_t i=0; i<l.variables.size(); ++i) {
                for(std::size_t j=i+1; j<l.variables.size(); ++j) {
                    const std::size_t var1 = l.variables[i].var;
                    const std::size_t var2 = l.variables[j].var;
                    adjacent_vars.insert({std::min(var1,var2), std::max(var1, var2)});
                }
            }
        }

        std::vector<std::size_t> adjacency_size(this->nr_variables(),0);
        for(const auto [i,j] : adjacent_vars) {
            ++adjacency_size[i];
            ++adjacency_size[j];
        }

        two_dim_variable_array<std::size_t> adjacency(adjacency_size.begin(), adjacency_size.end());
        std::fill(adjacency_size.begin(), adjacency_size.end(), 0);
        for(const auto e : adjacent_vars) {
            const auto [i,j] = e;
            assert(i<j);
            assert(adjacency_size[i] < adjacency[i].size());
            assert(adjacency_size[j] < adjacency[j].size());
            adjacency(i, adjacency_size[i]++) = j;
            adjacency(j, adjacency_size[j]++) = i;
        }

        return adjacency;
    }

    inline void ILP_input::reorder(const std::vector<std::size_t>& order)
    {
        assert(order.size() == this->nr_variables());
        std::vector<double> new_objective(this->nr_variables());
        for(std::size_t i=0; i<this->objective_.size(); ++i)
            if(order[i] < this->objective_.size())
                new_objective[i] = this->objective_[order[i]];
            else
                new_objective[i] = 0.0;
        std::swap(this->objective_, new_objective);

        std::vector<std::size_t> inverse_order(this->nr_variables());
        for(std::size_t i=0; i<order.size(); ++i)
            inverse_order[order[i]] = i;

        for(auto& [name, idx] : this->var_name_to_index_) {
            idx = inverse_order[idx];
        }

        std::vector<std::string> new_var_index_to_name(this->nr_variables());
        for(std::size_t i=0; i<this->var_index_to_name_.size(); ++i) {
            new_var_index_to_name[i] = this->var_index_to_name_[order[i]]; // TODO: use move operator
        }
        std::swap(new_var_index_to_name, this->var_index_to_name_);

        for(auto& l : this->linear_constraints_) {
            for(auto& x : l.variables) {
                x.var = inverse_order[x.var];
            }
            l.normalize(); 
        }
    }

    inline void ILP_input::reorder_bfs()
    {
        const auto adj = variable_adjacency_matrix();
        const auto order = bfs_ordering(adj);
        reorder(order);
    }

    inline void ILP_input::reorder_Cuthill_McKee()
    {
        const auto adj = variable_adjacency_matrix();
        const auto order = Cuthill_McKee(adj);
        reorder(order);
    }

}

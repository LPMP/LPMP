#pragma once

#include "config.hxx"
#include "bdd.h"
//#include "cuddObj.hh"
#include "hash_helper.hxx"
#include <iostream>
#include <unordered_map>
#include <numeric>
#include <tuple>

namespace LPMP {

    class bdd_converter {
        public:
            bdd_converter(BDD::bdd_mgr& bdd_mgr) : bdd_mgr_(bdd_mgr) 
        {}

            template<typename LEFT_HAND_SIDE_ITERATOR>
                BDD::node_ref convert_to_bdd(LEFT_HAND_SIDE_ITERATOR begin, LEFT_HAND_SIDE_ITERATOR end, const inequality_type ineq, const int right_hand_side);

            BDD::node_ref convert_to_bdd(const std::vector<int> coefficients, const inequality_type ineq, const int right_hand_side); 

        private:
            // returned vector has as its first element the right hand side, then follow the coefficients
            template<typename COEFF_ITERATOR>
                static std::tuple< std::vector<int>, inequality_type >
                normal_form(COEFF_ITERATOR begin, COEFF_ITERATOR end, const inequality_type ineq, const int right_hand_side);

            BDD::node_ref convert_to_bdd_impl(std::vector<int>& nf, const inequality_type ineq);

            BDD::bdd_mgr& bdd_mgr_;
            using constraint_cache_type = std::unordered_map<std::vector<int>,BDD::node_ref>;
            constraint_cache_type equality_cache;
            constraint_cache_type lower_equal_cache;
    };

    template<typename COEFF_ITERATOR>
        std::tuple< std::vector<int>, inequality_type >
        bdd_converter::normal_form(COEFF_ITERATOR begin, COEFF_ITERATOR end, const inequality_type ineq, const int right_hand_side)
        {
            assert(std::distance(begin,end) >= 1);
            int d = std::gcd(right_hand_side, *begin);
            for(auto it = begin+1; it != end; ++it)
                d = std::gcd(d, *it);

            std::vector<int> c;
            c.reserve(std::distance(begin, end) + 1);
            c.push_back(right_hand_side/d);
            for(auto it = begin; it != end; ++it)
                c.push_back(*it/d);

            if(ineq == inequality_type::greater_equal)
                for(auto& x : c)
                    x *= -1.0;

            return {c, ineq != inequality_type::greater_equal ? ineq : inequality_type::smaller_equal};
        }

    template<typename LEFT_HAND_SIDE_ITERATOR>
        BDD::node_ref bdd_converter::convert_to_bdd(LEFT_HAND_SIDE_ITERATOR begin, LEFT_HAND_SIDE_ITERATOR end, const inequality_type ineq, const int right_hand_side)
        {
            auto [nf, ineq_nf] = normal_form(begin, end, ineq, right_hand_side);
            return convert_to_bdd_impl(nf, ineq_nf); 
        }

    BDD::node_ref bdd_converter::convert_to_bdd(const std::vector<int> coefficients, const inequality_type ineq, const int right_hand_side)
        {
            return convert_to_bdd(coefficients.begin(), coefficients.end(), ineq, right_hand_side);
        }


    BDD::node_ref bdd_converter::convert_to_bdd_impl(std::vector<int>& nf, const inequality_type ineq)
        {
            assert(nf.size() > 0);
            const int right_hand_side = nf[0];
            if(nf.size() == 1) {
                switch(ineq) {
                    case inequality_type::equal: 
                        return right_hand_side == 0 ? bdd_mgr_.topsink() : bdd_mgr_.botsink();
                        //return right_hand_side == 0 ? bdd_mgr_.bddOne() : bdd_mgr_.bddZero();
                        break;
                    case inequality_type::smaller_equal:
                        return right_hand_side >= 0 ? bdd_mgr_.topsink() : bdd_mgr_.botsink();
                        break;
                    case inequality_type::greater_equal:
                        return right_hand_side <= 0 ? bdd_mgr_.topsink() : bdd_mgr_.botsink();
                        break;
                    default:
                        throw std::runtime_error("inequality type not supported");
                        break;
                }
            }

            // check, if constraint has already been seen. If so, retrieve
            switch(ineq) {
                case inequality_type::equal: 
                    {
                        auto cached = equality_cache.find(nf);
                        if(cached != equality_cache.end())
                            return cached->second;
                    }
                    break;
                case inequality_type::smaller_equal:
                    {
                        auto cached = lower_equal_cache.find(nf);
                        if(cached != lower_equal_cache.end())
                            return cached->second;
                    }
                    break;
                case inequality_type::greater_equal:
                    throw std::runtime_error("greater equal constraint not in normal form");
                    break;
                default:
                    throw std::runtime_error("inequality type not supported");
                    break;
            }

            // otherwise, recurse and build up constraint
            const int cur_coefficient = nf.back();
            nf.resize(nf.size()-1);

            BDD::node_ref bdd_0 = convert_to_bdd_impl(nf, ineq);
            // set first var to 1
            nf[0] -= cur_coefficient;
            BDD::node_ref bdd_1 = convert_to_bdd_impl(nf, ineq);
            nf[0] += cur_coefficient;

            BDD::node_ref cur_var = bdd_mgr_.projection(nf.size()-1);
            //auto bdd = cur_var.Ite(bdd_0, bdd_1);
            auto bdd = bdd_mgr_.ite_rec(cur_var, bdd_0, bdd_1);

            nf.push_back(cur_coefficient);
            // record bdd in cache
            switch(ineq) {
                case inequality_type::equal: 
                    equality_cache.insert(std::make_pair(nf,bdd));
                    break;
                case inequality_type::smaller_equal:
                    lower_equal_cache.insert(std::make_pair(nf,bdd));
                    break;
                case inequality_type::greater_equal:
                    throw std::runtime_error("greater equal constraint not in normal form");
                    break;
                default:
                    throw std::runtime_error("inequality type not supported");
                    break;
            } 

            return bdd;
        }

}

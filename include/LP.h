#pragma once

#include "config.hxx"
#include "factors_storage.hxx"
#include "messages_storage.hxx"
#include "message_passing_schedule.hxx"
#include "message_passing_weight_computation.hxx"
#include "lp_reparametrization.hxx"
#include "factor_container_interface.h"
#include <vector>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <functional>
#include <utility>
#include <limits>
#include <exception>
#include <unordered_map>
#include <unordered_set>
#include <atomic>
#include "template_utilities.hxx"
#include <assert.h>
#include "topological_sort.hxx"
#include <memory>
#include <iterator>
#include "two_dimensional_variable_array.hxx"
#include "union_find.hxx"
#include "serialization.hxx"
#include "tclap/CmdLine.h"

namespace LPMP {

template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename RECEIVE_MASK_ITERATOR>
void compute_pass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt, RECEIVE_MASK_ITERATOR receive_it)
{
   std::size_t n = std::distance(factorIt,factorItEnd);
   for(std::size_t i=0; i<n; ++i) {
      auto* f = *(factorIt + i);
      assert(f->FactorUpdated());
      f->UpdateFactor(*(omegaIt + i), *(receive_it + i));
   } 
}

template<typename FMC_TYPE>
class LP : public factors_storage<FMC_TYPE>, public messages_storage<FMC_TYPE> {
public:
   using FMC = FMC_TYPE;

   LP(TCLAP::CmdLine& cmd);
   LP(LP& o);

   void Begin(); // must be called after all messages and factors have been added
   void End() {};

   template<typename FACTOR_CONTAINER_TYPE, typename... ARGS>
   FACTOR_CONTAINER_TYPE* add_factor(ARGS&&... args);

   template<typename MESSAGE_CONTAINER_TYPE, typename LEFT_FACTOR, typename RIGHT_FACTOR, typename... ARGS>
   MESSAGE_CONTAINER_TYPE* add_message(LEFT_FACTOR* l, RIGHT_FACTOR* r, ARGS&&... args);

   //void ComputeWeights(const lp_reparametrization_mode m);
   void set_reparametrization(const lp_reparametrization r) { repam_mode_ = r; }
   lp_reparametrization get_repam_mode() const { return repam_mode_; }

   double LowerBound() const;
   void init_primal();
   double EvaluatePrimal() const;

   bool CheckPrimalConsistency() const;

   void ComputePass();
   void ComputeForwardPass();
   void ComputeBackwardPass();

   void ComputePassAndPrimal();
   void ComputeForwardPassAndPrimal();
   void ComputeBackwardPassAndPrimal();

   // compute pass with interleaved primal rounding on subset of potentials only. This can be used for horizon tracking and discrete tomography.
   template<typename FACTOR_ITERATOR, Direction DIRECTION>
   void ComputePassAndPrimal(FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end);

   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename RECEIVE_MASK_ITERATOR>
   void ComputePassAndPrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, OMEGA_ITERATOR omegaIt, RECEIVE_MASK_ITERATOR receive_mask_it);

   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename RECEIVE_MASK_ITERATOR>
   void ComputePass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt, RECEIVE_MASK_ITERATOR receive_it);

   struct message_passing_weight_storage 
   {
      //message_passing_weight_storage(std::initializer_list<> l) {assert(false);} // TODO: fill out
      weight_array omega_forward, omega_backward;
      receive_array receive_mask_forward, receive_mask_backward;
   };

   message_passing_weight_storage& get_message_passing_weight(const lp_reparametrization repam);

   double get_constant() const { return constant_; }

   void add_to_constant(const REAL x) { 
#pragma omp critical
      {
         constant_ += x; 
      }
   }

protected:

   //tsl::robin_map<lp_reparametrization, message_passing_weight_storage> message_passing_weights_;
   std::unordered_map<lp_reparametrization, message_passing_weight_storage> message_passing_weights_;
   lp_reparametrization repam_mode_ = lp_reparametrization(lp_reparametrization_mode::Undefined, 0.0);
   std::size_t rounding_iteration_ = 1;
   double constant_ = 0.0;

   TCLAP::ValueArg<std::string> reparametrization_type_arg_; // shared|residual|partition|overlapping_partition
   TCLAP::ValueArg<INDEX> inner_iteration_number_arg_;
   enum class reparametrization_type {shared,residual,partition,overlapping_partition};
   reparametrization_type reparametrization_type_ = reparametrization_type::shared;

   std::vector<bool> get_inconsistent_mask(const std::size_t no_fatten_rounds = 1);
   template<typename FACTOR_ITERATOR, typename FACTOR_MASK_ITERATOR>
   std::vector<FactorTypeAdapter*> get_masked_factors( FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end, FACTOR_MASK_ITERATOR factor_mask_begin, FACTOR_MASK_ITERATOR factor_mask_end);
   void reduce_optimization_factors();
};

template<typename FMC> 
LP<FMC>::LP(TCLAP::CmdLine& cmd)
: reparametrization_type_arg_("","reparametrizationType","message sending type: ", false, "shared", "{shared|residual|partition|overlapping_partition}", cmd)
, inner_iteration_number_arg_("","innerIteration","number of iterations in inner loop in partition reparamtrization, default = 5",false,5,&positiveIntegerConstraint,cmd)
{}

// make a deep copy of factors and messages. Adjust pointers to messages and factors
template<typename FMC>
LP<FMC>::LP(LP& o) // no const because of o.num_lp_threads_arg_.getValue() not being const!
  : reparametrization_type_arg_("","reparametrizationType","message sending type: ", false, o.reparametrization_type_arg_.getValue(), "{shared|residual|partition|overlapping_partition}" )
, inner_iteration_number_arg_("","innerIteration","number of iterations in inner loop in partition reparamtrization, default = 5",false,o.inner_iteration_number_arg_.getValue(),&positiveIntegerConstraint) 
{
  /*
  f_.reserve(o.f_.size());
  assert(false);
  std::unordered_map<FactorTypeAdapter*, FactorTypeAdapter*> factor_map; // translate addresses from o's factors to this' factors
  f_.reserve(o.f_.size());
  for(auto* f : o.f_) {
    auto* clone = f->clone();
    assert(false);
    //this->AddFactor(clone);
    factor_map.insert(std::make_pair(f, clone));
  }
  m_.reserve(o.m_.size());
  for(auto m : o.m_) {
    auto* left = m.left;
    auto* right = m.right;
    auto* left_clone = factor_map[left];
    auto* right_clone = factor_map[right];
    assert(false);
    //m_.push_back({left_clone, right_clone, m.sends_message_to_left, m.sends_message_to_right, m.receives_message_from_left, m.receives_message_from_right});
  }

  ordering_valid_ = o.ordering_valid_;

  message_passing_weights_ = o.message_passing_weights_;

  auto copy_factors = [&factor_map] (auto& factor_orig) {
     decltype(factor_orig) factor_copy;
     factor_copy.reserve(factor_orig.size());
     for(auto* f : factor_orig) {
        factor_copy.push_back( factor_map[f] );
     }
     return factor_copy;
  };

  forward_sorting = copy_factors(o.forward_sorting);
  backward_sorting = copy_factors(o.backward_sorting);

  forward_update_sorting = copy_factors(o.forward_update_sorting);
  backward_update_sorting = copy_factors(o.backward_update_sorting);

  auto copy_factor_rel = [&factor_map] (auto& factor_rel_orig) {
     decltype(factor_rel_orig) factor_rel_copy;
     factor_rel_copy.reserve(factor_rel_orig.size());
     for(auto f : factor_rel_orig) {
        factor_rel_copy.push_back( std::make_pair(factor_map[f.first], factor_map[f.second]) );
     }
     return factor_rel_orig; 
  };

  forward_pass_factor_rel_ = copy_factor_rel(o.forward_pass_factor_rel_);
  backward_pass_factor_rel_ = copy_factor_rel(o.backward_pass_factor_rel_);

  constant_ = o.constant_;
   */
}

template<typename FMC>
inline void LP<FMC>::Begin()
{
   repam_mode_ = lp_reparametrization();
   assert(repam_mode_.mode == lp_reparametrization_mode::Undefined);

   if(reparametrization_type_arg_.getValue() == "shared") {
     reparametrization_type_ = reparametrization_type::shared;
   } else if(reparametrization_type_arg_.getValue() == "residual") {
     reparametrization_type_ = reparametrization_type::residual;
   } else if(reparametrization_type_arg_.getValue() == "partition") {
     reparametrization_type_ = reparametrization_type::partition;
   } else if(reparametrization_type_arg_.getValue() == "overlapping_partition") {
     reparametrization_type_ = reparametrization_type::overlapping_partition;
   } else {
     throw std::runtime_error("reparamerization type not recognized");
   }
}

template<typename FMC>
template<typename FACTOR_CONTAINER_TYPE, typename... ARGS>
FACTOR_CONTAINER_TYPE* LP<FMC>::add_factor(ARGS&&... args)
{
   message_passing_weights_.clear();
   return factors_storage<FMC>::template add_factor<FACTOR_CONTAINER_TYPE>(std::forward<ARGS>(args)...);
}

template<typename FMC>
template<typename MESSAGE_CONTAINER_TYPE, typename LEFT_FACTOR, typename RIGHT_FACTOR, typename... ARGS>
MESSAGE_CONTAINER_TYPE* LP<FMC>::add_message(LEFT_FACTOR* l, RIGHT_FACTOR* r, ARGS&&... args)
{
   message_passing_weights_.clear();
   return messages_storage<FMC>::template add_message<MESSAGE_CONTAINER_TYPE>(l,r, std::forward<ARGS>(args)...);
}

template<typename FMC>
inline void LP<FMC>::ComputePass()
{
   ComputeForwardPass();
   ComputeBackwardPass();
}

template<typename FMC>
void LP<FMC>::ComputeForwardPass()
{
  auto mpw = get_message_passing_weight(repam_mode_);
  auto [forward_sorting, forward_update_sorting] = this->get_sorted_factors(Direction::forward);
  ComputePass(forward_update_sorting.begin(), forward_update_sorting.end(), mpw.omega_forward.begin(), mpw.receive_mask_forward.begin()); 
}

template<typename FMC>
void LP<FMC>::ComputeBackwardPass()
{
  auto mpw = get_message_passing_weight(repam_mode_);
  auto [backward_sorting, backward_update_sorting] = this->get_sorted_factors(Direction::backward);
  ComputePass(backward_update_sorting.begin(), backward_update_sorting.end(), mpw.omega_backward.begin(), mpw.receive_mask_backward.begin()); 
}

template<typename FMC>
void LP<FMC>::ComputeForwardPassAndPrimal()
{
  auto mpw = get_message_passing_weight(repam_mode_);
  auto [forward_sorting, forward_update_sorting] = this->get_sorted_factors(Direction::forward);
  ComputePassAndPrimal(forward_update_sorting.begin(), forward_update_sorting.end(), mpw.omega_forward.begin(), mpw.receive_mask_forward.begin());
}

template<typename FMC>
void LP<FMC>::ComputeBackwardPassAndPrimal()
{
  auto mpw = get_message_passing_weight(repam_mode_);
  auto [backward_sorting, backward_update_sorting] = this->get_sorted_factors(Direction::backward);
  ComputePassAndPrimal(backward_update_sorting.begin(), backward_update_sorting.end(), mpw.omega_backward.begin(), mpw.receive_mask_backward.begin()); 
}

template<typename FMC>
void LP<FMC>::ComputePassAndPrimal()
{
  ComputeForwardPassAndPrimal();
  ComputeBackwardPassAndPrimal();
}

template<typename FMC>
template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename RECEIVE_MASK_ITERATOR>
void LP<FMC>::ComputePass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt, RECEIVE_MASK_ITERATOR receive_it)
{
    const std::size_t n = std::distance(factorIt, factorItEnd);
    //#pragma omp parallel for schedule(static)
    if(reparametrization_type_ == reparametrization_type::shared || reparametrization_type_ == reparametrization_type::partition || reparametrization_type_ == reparametrization_type::overlapping_partition) {
        for(std::size_t i=0; i<n; ++i) {
            auto* f = *(factorIt + i);
            assert(f->FactorUpdated());
            f->UpdateFactor(*(omegaIt + i), *(receive_it + i));
        }
    } else if(reparametrization_type_ == reparametrization_type::residual) {
        for(std::size_t i=0; i<n; ++i) {
            auto* f = *(factorIt + i);
            f->update_factor_residual(*(omegaIt + i), *(receive_it + i));
        }
    } else {
       throw std::runtime_error("reparametrization type not recognized");
    }
}

template<typename FMC>
typename LP<FMC>::message_passing_weight_storage& 
LP<FMC>::get_message_passing_weight(const lp_reparametrization repam)
{
   auto it = message_passing_weights_.find(repam);
   if(it != message_passing_weights_.end()) {
      return it->second;
   } else {
      auto [forward_sorting, forward_update_sorting] = this->get_sorted_factors(Direction::forward);
      auto [backward_sorting, backward_update_sorting] = this->get_sorted_factors(Direction::backward);
      const auto leave_percentage = repam.leave_percentage;
      if(repam.mode == lp_reparametrization_mode::Anisotropic) {
         auto [omega_forward, receive_mask_forward] = compute_anisotropic_weights(forward_sorting.begin(), forward_sorting.end(), leave_percentage);
         auto [omega_backward, receive_mask_backward] = compute_anisotropic_weights(backward_sorting.begin(), backward_sorting.end(), leave_percentage);
         message_passing_weight_storage mpw {omega_forward, omega_backward, receive_mask_forward, receive_mask_backward};
         message_passing_weights_.insert( std::make_pair(repam, mpw) ); // TODO: move 
      } else if(repam.mode == lp_reparametrization_mode::Uniform) {
         auto omega_forward = compute_isotropic_weights(forward_sorting.begin(), forward_sorting.end(), leave_percentage);
         omega_valid(forward_update_sorting.begin(), forward_update_sorting.end(), omega_forward);
         auto omega_backward = compute_isotropic_weights(backward_sorting.begin(), backward_sorting.end(), leave_percentage);
         omega_valid(backward_update_sorting.begin(), backward_update_sorting.end(), omega_backward);
         auto receive_mask_forward = compute_full_receive_mask(forward_sorting.begin(), forward_sorting.end());
         receive_mask_valid(forward_update_sorting.begin(), forward_update_sorting.end(), receive_mask_forward);
         auto receive_mask_backward = compute_full_receive_mask(backward_sorting.begin(), backward_sorting.end());
         receive_mask_valid(backward_update_sorting.begin(), backward_update_sorting.end(), receive_mask_backward);
         message_passing_weight_storage mpw {omega_forward, omega_backward, receive_mask_forward, receive_mask_backward};
         message_passing_weights_.insert( std::make_pair(repam, mpw) ); // TODO: move 
      } else if(repam.mode == lp_reparametrization_mode::Anisotropic2) {
         auto [omega_forward, receive_mask_forward] = compute_anisotropic_weights_2(forward_sorting.begin(), forward_sorting.end(), leave_percentage);
         auto [omega_backward, receive_mask_backward] = compute_anisotropic_weights_2(backward_sorting.begin(), backward_sorting.end(), leave_percentage);
         message_passing_weight_storage mpw {omega_forward, omega_backward, receive_mask_forward, receive_mask_backward};
         message_passing_weights_.insert( std::make_pair(repam, mpw) ); // TODO: move 
      } else {
         throw std::runtime_error("no reparametrization mode set");
      }
      return message_passing_weights_[repam];
   } 
}

// check whether messages constraints are satisfied
template<typename FMC>
bool LP<FMC>::CheckPrimalConsistency() const
{
   std::atomic<bool> consistent=true; // or use std::atomic<bool> for multithreading?

   this->for_each_factor([&consistent](auto& f) {
         if(consistent && !f.check_primal_consistency()) {
            consistent = false;
         }
   });

   if(debug()) { std::cout << "primal solution consistent: " << (consistent ? "true" : "false") << "\n"; }
   return consistent;
}

template<typename FMC>
double LP<FMC>::LowerBound() const
{
    double lb = constant_;
    this->for_each_factor([&lb](const auto& f) { 
       const double delta = f.LowerBound();
       assert(delta < std::numeric_limits<double>::max());
#pragma omp critical
       {
         lb += delta;
       }
       assert(std::isfinite(lb));
    });

    //assert(std::abs(lb - std::accumulate(factors_storage<FMC>::cbegin(), factors_storage<FMC>::cend(), [](const double lb, const auto* f) { return lb + f->LowerBound; }, 0.0)) <= eps);
    return lb;
}

template<typename FMC>
void LP<FMC>::init_primal()
{
    this->for_each_factor([](auto& f) { 
       f.init_primal();
    }); 
}

template<typename FMC>
double LP<FMC>::EvaluatePrimal() const
{
    const bool consistent = CheckPrimalConsistency();
    if(consistent == false)
       return std::numeric_limits<REAL>::infinity();

    double cost = constant_;
    this->for_each_factor([&cost](const auto& f) {
          const double delta = f.EvaluatePrimal();
          //assert(delta < std::numeric_limits<double>::max()); // need not hold true, since we can have infeasible solutions
#pragma omp critical
          {
            cost += delta;
          }
    });
    if(debug())
       std::cout << "primal cost = " << cost << "\n";
    return cost;
}


template<typename FMC>
template<typename FACTOR_ITERATOR, Direction DIRECTION>
void LP<FMC>::ComputePassAndPrimal(FACTOR_ITERATOR factor_begin, const FACTOR_ITERATOR factor_end)
{
    // filter out factors that that are not needed
    std::unordered_set<FactorTypeAdapter*> factor_set;
    for(auto f_it=factor_begin; f_it!=factor_end; ++f_it) { factor_set.insert(*f_it); }
    assert(factor_set.size() == std::distance(factor_begin, factor_end));
    auto filter_predicate = [&factor_set](FactorTypeAdapter* f) { return factor_set.count(f) == 0; };

    // forward
    std::vector<FactorTypeAdapter*> filtered_factors;
    std::vector<FactorTypeAdapter*> filtered_factors_update;
    if(DIRECTION == Direction::forward) {
        auto [forward_sorting, forward_update_sorting] = this->get_sorted_factors(Direction::forward);
        std::remove_copy_if(forward_sorting.begin(), forward_sorting.end(), std::back_inserter(filtered_factors), filter_predicate);
        std::remove_copy_if(forward_update_sorting.begin(), forward_update_sorting.end(), std::back_inserter(filtered_factors_update), filter_predicate);
    } else {
        assert(DIRECTION == Direction::backward);
        auto [backward_sorting, backward_update_sorting] = this->get_sorted_factors(Direction::backward);
        std::remove_copy_if(backward_update_sorting.begin(), backward_update_sorting.end(), std::back_inserter(filtered_factors_update), filter_predicate);
        std::remove_copy_if(backward_sorting.begin(), backward_sorting.end(), std::back_inserter(filtered_factors), filter_predicate);
    }

    assert(std::distance(factor_begin, factor_end) == filtered_factors.size());

    auto [omega, receive_mask] = compute_anisotropic_weights( filtered_factors.begin(), filtered_factors.end(), 0.0);

    assert(filtered_factors_update.size() == omega.size());
    assert(filtered_factors_update.size() == receive_mask.size());

    // prune masks: set weights going to factors outside to zero
    for(std::size_t i=0; i<filtered_factors_update.size(); ++i) {
        auto messages = filtered_factors_update[i]->get_messages();
        std::size_t j_omega = 0;
        std::size_t j_receive = 0;
        for(std::size_t j=0; j<messages.size(); ++j) {
            auto* adjacent_factor = messages[j].adjacent_factor;
            const bool adjacent_factor_in = !filter_predicate(adjacent_factor);
            if(message_passing_schedule_factor_view::sends_message_to_adjacent_factor(messages[j].mps) && !adjacent_factor_in) {
                receive_mask[i][j_receive++] = 0;
            }
            if(message_passing_schedule_factor_view::receives_message_from_adjacent_factor(messages[j].mps) && !adjacent_factor_in) {
                omega[i][j_omega++] = 0.0;
            }
        } 
    }

    ComputePassAndPrimal(filtered_factors_update.begin(), filtered_factors_update.end(), omega.begin(), receive_mask.begin());
}

template<typename FMC>
template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename RECEIVE_MASK_ITERATOR>
void LP<FMC>::ComputePassAndPrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, OMEGA_ITERATOR omegaIt, RECEIVE_MASK_ITERATOR receive_mask_it)
{
   rounding_iteration_++;

   for(std::size_t i=0; i<std::distance(factorIt, factorEndIt); ++i) {
      auto* f = *(factorIt+i);
      assert(f->FactorUpdated());
      f->UpdateFactorPrimal(*(omegaIt + i), *(receive_mask_it + i), rounding_iteration_);
   }
}

template<typename FMC>
std::vector<bool> LP<FMC>::get_inconsistent_mask(const std::size_t no_fatten_rounds)
{
  std::vector<bool> inconsistent_mask;
  inconsistent_mask.reserve(this->number_of_factors());
  
  // check for locally non-optimal factors and check for violated messages
  for(auto factor_it = factors_storage<FMC>::begin(); factor_it != factors_storage<FMC>::end(); ++factor_it) {
    auto* f = *factor_it;
    assert(f->EvaluatePrimal() < std::numeric_limits<REAL>::infinity());
    if(f->LowerBound() < f->EvaluatePrimal() - eps || !f->check_primal_consistency()) {
      inconsistent_mask.push_back(true);
    } else {
      inconsistent_mask.push_back(false);
    } 
  }

  // fatten the region
  auto fatten = [&]() {
     for(auto msg_it = messages_storage<FMC>::begin(); msg_it != messages_storage<FMC>::end(); ++msg_it) {
        auto* l = msg_it->left;
        auto l_index = this->get_factor_index(l);
        auto* r = msg_it->right;
        auto r_index = this->get_factor_index(r);

        if(inconsistent_mask[l_index] == true || inconsistent_mask[r_index] == true) {
           inconsistent_mask[l_index] = true;
           inconsistent_mask[r_index] = true;
        }
    }
  };

  for(std::size_t iter=0; iter<no_fatten_rounds; ++iter) {
    fatten();
  }

  if(debug()) {
    std::cout << "\% inconsistent factors = " << std::count(inconsistent_mask.begin(), inconsistent_mask.end(), true)/REAL(this->number_of_factors()) << "\n";
  }

  return inconsistent_mask;
}

template<typename FMC>
template<typename FACTOR_ITERATOR, typename FACTOR_MASK_ITERATOR>
std::vector<FactorTypeAdapter*> LP<FMC>::get_masked_factors(
    FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end,
    FACTOR_MASK_ITERATOR factor_mask_begin, FACTOR_MASK_ITERATOR factor_mask_end)
{
  assert(std::distance(factor_mask_begin, factor_mask_end) == this->number_of_factors());
  
  std::vector<FactorTypeAdapter*> factors;
  for(auto f_it=factor_begin; f_it!=factor_end; ++f_it) {
    const auto f_index = this->get_factor_address(*f_it);
    if(factor_mask_begin[f_index]) {
      factors.push_back(*f_it);
    }
  }

  return factors;
}

} // end namespace LPMP

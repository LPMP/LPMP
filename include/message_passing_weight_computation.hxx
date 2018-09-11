#ifndef LPMP_MESSAGE_PASSING_WEIGHT_COMPUTATION_HXX
#define LPMP_MESSAGE_PASSING_WEIGHT_COMPUTATION_HXX

#include "two_dimensional_variable_array.hxx"
#include "message_passing_schedule.hxx"
#include "factor_container_interface.h"
#include <tuple>
#include <vector>
#include <unordered_map>
#include <cassert>

namespace LPMP {

using weight_array = two_dim_variable_array<double>;
using weight_slice = two_dim_variable_array<double>::ArrayAccessObject;
using receive_array = two_dim_variable_array<unsigned char>;
using receive_slice = two_dim_variable_array<unsigned char>::ArrayAccessObject;

template<typename ITERATOR>
weight_array allocate_omega(ITERATOR factor_begin, ITERATOR factor_end)
{
   std::vector<std::size_t> omega_size;
   omega_size.reserve(std::distance(factor_begin, factor_end));

   for(auto factor_it=factor_begin; factor_it != factor_end; ++factor_it) {
       if((*factor_it)->FactorUpdated()) {
           omega_size.push_back((*factor_it)->no_send_messages());
       }
   }

   weight_array omega(omega_size.begin(), omega_size.end());

   assert(omega.size() == omega_size.size());
   for(std::size_t i=0; i<omega.size(); ++i) { assert(omega[i].size() == omega_size[i]); }

   return omega;
}

template<typename ITERATOR>
receive_array allocate_receive_mask(ITERATOR factor_begin, ITERATOR factor_end)
{
   std::vector<std::size_t> receive_mask_size;
   receive_mask_size.reserve(std::distance(factor_begin, factor_end));

   for(auto factor_it=factor_begin; factor_it != factor_end; ++factor_it) {
       if((*factor_it)->FactorUpdated()) {
           receive_mask_size.push_back((*factor_it)->no_receive_messages());
       }
   }

   receive_array receive_mask(receive_mask_size.begin(), receive_mask_size.end());

   assert(receive_mask.size() == receive_mask_size.size());
   for(std::size_t i=0; i<receive_mask.size(); ++i) { assert(receive_mask[i].size() == receive_mask_size[i]); }

   return receive_mask;
}

// compute isotropic weights so as to help decoding for obtaining primal solutions and tightening
// leave_percentage signals how much weight to leave in sending factor. Important for rounding and tightening.
// TODO: possibly pass update factor list (i.e. only those factors which are updated)
template<typename FACTOR_ITERATOR>
weight_array compute_isotropic_weights
(
 FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end,
 const double leave_percentage)
{
   assert(leave_percentage >= 0.0 && leave_percentage < 1.0);

   auto omega = allocate_omega(factor_begin, factor_end);

   std::size_t c=0;
   for(auto factor_it=factor_begin; factor_it != factor_end; ++factor_it) {
     if((*factor_it)->FactorUpdated()) {
       std::size_t k=0;
       auto msgs = (*factor_it)->get_messages();
       const auto weight = (1.0-leave_percentage)/double((*factor_it)->no_send_messages());
       for(auto msg : msgs) {
           if(message_passing_schedule_factor_view::sends_message_to_adjacent_factor(msg.mps)) {
               omega[c][k] = weight;
               ++k;
         }
       }
       assert(k == omega[c].size());
       c++;
     }
   }
   assert(c == omega.size());

   return omega;
}

template<typename FACTOR_ITERATOR>
std::unordered_map<FactorTypeAdapter*, std::size_t> get_factor_indices(FACTOR_ITERATOR f_begin, FACTOR_ITERATOR f_end)
{
    std::unordered_map<FactorTypeAdapter*, std::size_t> indices;
    indices.reserve( std::distance(f_begin, f_end) );
    std::size_t i=0;
    for(auto f_it=f_begin; f_it!=f_end; ++f_it, ++i) {
        indices.insert( {*f_it, i} );
    }

    return indices;
}

template<typename FACTOR_ITERATOR>
std::tuple<weight_array,receive_array> 
compute_anisotropic_weights(FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end, const double leave_percentage)
{
   assert(leave_percentage >= 0.0 && leave_percentage < 1.0);

   const auto n = std::distance(factor_begin,factor_end);

   auto factor_address_to_sorted_index = get_factor_indices(factor_begin, factor_end); 

   // compute the following numbers: 
   // 1) #{factors after current one, to which messages are sent from current factor}
   // 2) #{factors after current one, which receive messages from current one}
   std::vector<std::size_t> no_later_factors_receiving_message_from(n, 0); // number of factors later than current one receiving message from current one
   std::vector<std::size_t> last_factor_receiving_message_from(n, 0); // what is the last (in the order given by factor iterator) factor that receives a message?
   std::vector<std::size_t> first_factor_receiving_message_from(n, std::numeric_limits<std::size_t>::max()); // what is the last (in the order given by factor iterator) factor that receives a message?
   std::vector<std::size_t> first_factor_sending_message_to(n, std::numeric_limits<std::size_t>::max());

   bool connected_factor_not_in_list = false;

   for(auto f_it=factor_begin; f_it!=factor_end; ++f_it) {
      assert(factor_address_to_sorted_index.count(*f_it) > 0);
      const auto f_index = factor_address_to_sorted_index[*f_it];
      const auto messages = (*f_it)->get_messages();
      for(const auto m : messages) {
         if(factor_address_to_sorted_index.count(m.adjacent_factor) > 0) {
            const auto adjacent_index = factor_address_to_sorted_index[m.adjacent_factor]; 
            if(message_passing_schedule_factor_view::receives_message_to_adjacent_factor(m.mps) ) {
               last_factor_receiving_message_from[f_index] = std::max(last_factor_receiving_message_from[f_index], adjacent_index);
               first_factor_receiving_message_from[f_index] = std::min(first_factor_receiving_message_from[f_index], adjacent_index);
               if(adjacent_index > f_index) {
                  no_later_factors_receiving_message_from[f_index]++;
               }
            }
            if(message_passing_schedule_factor_view::sends_message_to_adjacent_factor(m.mps)) {
               first_factor_sending_message_to[adjacent_index] = std::min(f_index, first_factor_sending_message_to[adjacent_index]);
            } else {
               connected_factor_not_in_list = true;
            }
         }
      }
   }

   // now take into account factors that are not iterated over, but from which a factor that is iterated over may send and another can receive.
   // It still makes sense to send and receive from such factors
   std::unordered_map<FactorTypeAdapter*, std::size_t> min_adjacent_sending;
   min_adjacent_sending.reserve(n);
   std::unordered_map<FactorTypeAdapter*, std::size_t> max_adjacent_receiving;
   max_adjacent_receiving.reserve(n);


   // there are factors connected to the ones given to which we also send messages. We do so only if we can afterwards receive messages from those later.
   if(connected_factor_not_in_list) {

       // get vector of factors that are (i) not in iteration list and (ii) are connected to two or more factors in iteration list.
       std::unordered_map<FactorTypeAdapter*, std::size_t> adjacent_factors;
       for(auto f_it=factor_begin; f_it!=factor_end; ++f_it) {
           for(auto* f : (*f_it)->get_adjacent_factors()) {
               if(factor_address_to_sorted_index.count(f) == 0) {
                   // is adjacent_factors initialized with value 0?
                   adjacent_factors[f]++;
               }
           }
       }

       for(auto e : adjacent_factors) {
           if(e.second >= 2) {
               FactorTypeAdapter* f = e.first;
               auto& min_adjacent_sending_index = min_adjacent_sending[f];
               auto& max_adjacent_receiving_index = max_adjacent_receiving[f] = 0;
               min_adjacent_sending_index = std::numeric_limits<std::size_t>::max();
               max_adjacent_receiving_index = 0;
               const auto adjacent_factors = e.first->get_messages();
               for(const auto f : adjacent_factors) {
                   if(factor_address_to_sorted_index.count(f.adjacent_factor)) {
                       const auto adjacent_index = factor_address_to_sorted_index[f.adjacent_factor];
                       if(message_passing_schedule_factor_view::sends_message_from_adjacent_factor(f.mps)) {
                           min_adjacent_sending_index = std::min(adjacent_index, min_adjacent_sending_index);
                       }
                       if(message_passing_schedule_factor_view::receives_message_to_adjacent_factor(f.mps)) {
                           max_adjacent_receiving_index = std::max(adjacent_index, max_adjacent_receiving_index);
                       } 
                   }
               }
           }
       } 
   }

   auto omega = allocate_omega(factor_begin, factor_end);
   auto receive_mask = allocate_receive_mask(factor_begin, factor_end);

   auto receives_msg = [&](FactorTypeAdapter* factor, const std::size_t factor_index, auto m) -> bool
   { 
       auto* adjacent_factor = m.adjacent_factor;
       assert(adjacent_factor != factor);
       assert(factor_address_to_sorted_index.count(factor) > 0 && factor_address_to_sorted_index[factor] == factor_index);
       assert(message_passing_schedule_factor_view::receives_message_from_adjacent_factor(m.mps) == true);
       if(factor_address_to_sorted_index.count(adjacent_factor) > 0) {
           const auto adjacent_factor_index = factor_address_to_sorted_index[adjacent_factor];
           if(adjacent_factor_index < factor_index && adjacent_factor->FactorUpdated())  { return true; }
           if(first_factor_receiving_message_from[adjacent_factor_index] < factor_index) { return true; }
           return false;
       } else {
           assert(connected_factor_not_in_list);
           const auto min_adjacent_sending_index = min_adjacent_sending[adjacent_factor];
           if(min_adjacent_sending_index < factor_index) { return true; }
           return false;
       } 
   };

   auto sends_msg = [&](FactorTypeAdapter* factor, const std::size_t factor_index, auto m) -> bool
   { 
       auto* adjacent_factor = m.adjacent_factor;
       assert(adjacent_factor != factor);
       assert(factor_address_to_sorted_index.count(factor) > 0 && factor_address_to_sorted_index[factor] == factor_index);
       assert(message_passing_schedule_factor_view::sends_message_to_adjacent_factor(m.mps) == true);
       if(factor_address_to_sorted_index.count(adjacent_factor) > 0) {
           const auto adjacent_factor_index = factor_address_to_sorted_index[adjacent_factor];
           //if(m.adjacent_factor_receives && receives_msg(adjacent_factor, adjacent_factor_index, m.reverse(factor))) { return false; }
           if(factor_index < adjacent_factor_index && adjacent_factor->FactorUpdated()) { return true; }
           if(last_factor_receiving_message_from[adjacent_factor_index] > factor_index) { return true; }
           return false;
       } else {
           assert(connected_factor_not_in_list);
           const auto max_adjacent_receiving_index = max_adjacent_receiving[adjacent_factor];
           if(factor_index < max_adjacent_receiving_index) { return true; }
           return false; 
       }
   };


   {
      std::size_t c=0;
      for(auto f_it=factor_begin; f_it!=factor_end; ++f_it) {
          auto* factor = *f_it;
          if(factor->FactorUpdated()) {
              const auto factor_index = factor_address_to_sorted_index[factor];
              std::size_t k_send = 0;
              std::size_t k_receive = 0;
               
              // indicate which messages are sent and received
              const auto msgs = factor->get_messages();
              for(auto m : msgs) {
                  if(message_passing_schedule_factor_view::sends_message_to_adjacent_factor(m.mps)) {
                      if(sends_msg(factor, factor_index, m)) {
                          omega[c][k_send] = 1.0;
                      } else {
                          omega[c][k_send] = 0.0;
                      } 
                      ++k_send; 
                  }

                  if(message_passing_schedule_factor_view::receives_message_from_adjacent_factor(m.mps)) {
                      if(receives_msg(factor, factor_index, m)) {
                          receive_mask[c][k_receive] = 1; 
                      } else {
                          receive_mask[c][k_receive] = 0; 
                      }
                      ++k_receive; 
                  } 
              }

              assert(k_receive == factor->no_receive_messages());
              assert(k_send == factor->no_send_messages());

              // set omega reweighting values
              const std::size_t no_send_messages_anisotropic = std::count(omega[c].begin(), omega[c].end(), 1.0);
              const auto no_send_messages = factor->no_send_messages();
              const auto leave_weight = [&]() {
                  if(no_later_factors_receiving_message_from[factor_index] > 0) return 1.0;
                  if(no_send_messages - no_send_messages_anisotropic  > no_send_messages_anisotropic) return 1.0;
                  return 0.0;
              }();

              // higher weight than SRMP chooses, that still leaves something in the factor, exactly when SRMP does so.
              //const double weight = (1.0-leave_percentage) / double(leave_weight + no_send_messages_anisotropic);

              // srmp option:
              const double srmp_weight = (1.0-leave_percentage)/double(no_later_factors_receiving_message_from[factor_index] + std::max(no_send_messages_anisotropic, no_send_messages - no_send_messages_anisotropic));

              if(no_send_messages_anisotropic > 0) {
                  for(auto& x : omega[c]) { if(x > 0) { x *= srmp_weight; } }
              }

              assert(std::accumulate(omega[c].begin(), omega[c].end(), 0.0) <= 1.0 + eps);
              ++c;
          }
      }
      assert(c == omega.size());
   }

   // TODO: perform first check in LP class
   // check whether all messages were added to m_. Possibly, this can be automated: Traverse all factors, get all messages, add them to m_ and avoid duplicates along the way.
   //assert(2*m_.size() == std::accumulate(f_.begin(), f_.end(), 0, [](auto sum, auto* f){ return sum + f->no_messages(); }));
   for(std::size_t i=0; i<omega.size(); ++i) {
      assert(std::accumulate(omega[i].begin(), omega[i].end(), 0.0) <= 1.0 + eps);
   }

   return {omega,receive_mask};
}

template<typename FACTOR_ITERATOR>
std::tuple<weight_array, receive_array> 
compute_anisotropic_weights_2(FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end, const double leave_percentage)
{
   auto factor_address_to_sorted_index = get_factor_indices(factor_begin, factor_end); 

   auto omega = allocate_omega(factor_begin, factor_end);
   auto receive_mask = allocate_receive_mask(factor_begin, factor_end); 

   std::vector<std::size_t> no_send_messages_later(std::distance(factor_begin,factor_end), 0);
   /*
   for(std::size_t i=0; i<m_.size(); ++i) {
      auto* f_left = m_[i].left;
      const std::size_t f_index_left = factor_address_to_index_[f_left];
      const std::size_t index_left = f_sorted_inverse[f_index_left];
      auto* f_right = m_[i].right;
      const std::size_t f_index_right = factor_address_to_index_[f_right];
      const std::size_t index_right = f_sorted_inverse[f_index_right];
      
      if(m_[i].sends_message_to_right && index_left < index_right) {
        no_send_messages_later[index_left]++;
      }
      if(m_[i].sends_message_to_left && index_right < index_left) {
        no_send_messages_later[index_right]++;
      }
   }

   auto omega = allocate_omega(factor_begin, factor_end);
   auto receive_mask = allocate_receive_mask(factor_begin, factor_end);

   {
      std::size_t c=0;
      for(auto it=factor_begin; it!=factor_end; ++it) {
         const std::size_t i = std::distance(factor_begin, it);
         assert(i == f_sorted_inverse[ factor_address_to_index_[*it] ]);
         if((*it)->FactorUpdated()) {
            std::size_t k_send=0;
            std::size_t k_receive=0;
            auto msgs = (*it)->get_messages();
            for(auto msg_it : msgs) {
                auto* f_connected = msg_it.adjacent_factor;
                const std::size_t j = f_sorted_inverse[ factor_address_to_index_[f_connected] ];
                if(msg_it.sends_to_adjacent_factor) {
                    assert(i != j);
                    if(i<j) {
                        omega[c][k_send] = (1.0-leave_percentage)/double(no_send_messages_later[i]);
                    } else {
                        omega[c][k_send] = 0.0;
                    } 
                    ++k_send;
                }
                if(msg_it.receives_from_adjacent_factor) {
                    if(j<i) {
                        receive_mask[c][k_receive] = 1;
                    } else {
                        receive_mask[c][k_receive] = 0;
                    }
                    ++k_receive;
               }
            }
            ++c;
         }
      }
   }
   */

   omega_valid(factor_begin, factor_end, omega);
   return {omega, receive_mask};
}

template<typename FACTOR_ITERATOR>
receive_array compute_full_receive_mask(FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end)
{
    std::vector<std::size_t> mask_size;
    mask_size.reserve(std::distance(factor_begin, factor_end));
    for(auto it=factor_begin; it!=factor_end; ++it) {
        if((*it)->FactorUpdated()) {
            mask_size.push_back( (*it)->no_receive_messages() );
        } 
    }
    receive_array receive_mask(mask_size.begin(), mask_size.end());
    for(std::size_t i=0; i<receive_mask.size(); ++i) {
        for(std::size_t j=0; j<receive_mask[i].size(); ++j) {
            receive_mask(i,j) = true;
        }
    }

    return receive_mask;
}

template<typename FACTOR_ITERATOR>
void omega_valid(FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end, const weight_array& omega)
{
   auto factor_it = factor_begin;
   std::size_t i=0;
   for(; factor_it!=factor_end; ++factor_it) {
      if((*factor_it)->FactorUpdated()) {
         assert(*std::min_element(omega[i].begin(), omega[i].end()) >= 0.0);
         assert(std::accumulate(omega[i].begin(), omega[i].end(), 0.0) <= 1.0 + eps);
         assert(omega[i].size() == factor_begin[i]->no_send_messages());
         ++i;
      }
   }
   assert(i == omega.size());
   assert(factor_it == factor_end);
}

template<typename FACTOR_ITERATOR>
void receive_mask_valid(FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end, const receive_array& receive_mask)
{
   auto factor_it = factor_begin;
   std::size_t i=0;
   for(; factor_it!=factor_end; ++factor_it) {
      if((*factor_it)->FactorUpdated()) {
         assert(*std::min_element(receive_mask[i].begin(), receive_mask[i].end()) >= 0);
         assert(*std::max_element(receive_mask[i].begin(), receive_mask[i].end()) <= 1);
         assert(receive_mask[i].size() == factor_begin[i]->no_receive_messages());
         ++i;
      }
   }
   assert(i == receive_mask.size());
   assert(factor_it == factor_end);
}

} // namespace LPMP

#endif // LPMP_MESSAGE_PASSING_WEIGHT_COMPUTATION_HXX

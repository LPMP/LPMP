#pragma once

#include <vector>
#include <array>
#include <unordered_map>
#include <tuple>
#include <cassert>
#include <tsl/robin_map.h>
#include "meta/meta.hpp"
#include "topological_sort.hxx"
#include "factor_container_interface.h"

namespace LPMP {

template<typename FMC>
class factors_storage 
{
public:
   ~factors_storage() { for(auto* f : factors_) { delete f; } }

   std::size_t number_of_factors() const { return factors_.size(); }

   FactorTypeAdapter* get_factor(const std::size_t i) const { assert(i<number_of_factors()); return factors_[i]; }

   template<typename FACTOR_CONTAINER_TYPE, typename... ARGS>
   FACTOR_CONTAINER_TYPE* add_factor(ARGS&&... args);

   void add_factor_relation(FactorTypeAdapter* f1, FactorTypeAdapter* f2); // indicate that factor f1 comes before factor f2
   void add_forward_pass_factor_relation(FactorTypeAdapter* f1, FactorTypeAdapter* f2);
   void add_backward_pass_factor_relation(FactorTypeAdapter* f1, FactorTypeAdapter* f2);

   std::tuple<std::vector<FactorTypeAdapter*>&, std::vector<FactorTypeAdapter*>&> get_sorted_factors(const Direction d);

   template<typename FUNC> void for_each_factor(FUNC f);
   template<typename FUNC> void for_each_factor(FUNC f) const;

   auto begin() { return factors_.begin(); }
   auto end() { return factors_.end(); }

   auto cbegin() const { return factors_.cbegin(); }
   auto cend() const { return factors_.cend(); }

   std::size_t get_factor_index(const FactorTypeAdapter* f) const;

private:
   template<typename FACTOR_CONTAINER_TYPE>
   static constexpr std::size_t factors_tuple_index();

   // (ordering, ordering of factors that are updated)
   std::tuple<std::vector<FactorTypeAdapter*>, std::vector<FactorTypeAdapter*>>
      sort_factors(const std::vector<std::array<FactorTypeAdapter*,2>>& factor_rel);

   std::vector<FactorTypeAdapter*> factors_;
   tsl::robin_map<const FactorTypeAdapter*,std::size_t> factor_address_to_index_;
   //std::unordered_map<const FactorTypeAdapter*,std::size_t> factor_address_to_index_;
   std::vector<std::array<FactorTypeAdapter*,2>> forward_pass_factor_relation_, backward_pass_factor_relation_;
   std::vector<FactorTypeAdapter*> forward_pass_factor_ordering_, backward_pass_factor_ordering_;
   std::vector<FactorTypeAdapter*> forward_pass_factor_update_ordering_, backward_pass_factor_update_ordering_;

   struct vector_of_pointers {
      template<class T> 
         using invoke = typename std::vector<T*>;
   };
   using factors_vector_list = meta::transform< typename FMC::FactorList, vector_of_pointers >;
   using factors_tuple_type = meta::apply<meta::quote<std::tuple>, factors_vector_list>;

   factors_tuple_type factors_tuple_;
};

template<typename FMC>
template<typename FACTOR_CONTAINER_TYPE, typename... ARGS>
FACTOR_CONTAINER_TYPE* factors_storage<FMC>::add_factor(ARGS&&... args)
{ 
   auto* f = new FACTOR_CONTAINER_TYPE(std::forward<ARGS>(args)...);
   assert(factor_address_to_index_.size() == factors_.size());
   factors_.push_back(f);

   assert(factor_address_to_index_.count(f) == 0);
   factor_address_to_index_.insert(std::make_pair(f,factors_.size()-1));
   assert(factor_address_to_index_.find(f)->second == factors_.size()-1);

   constexpr auto factor_idx = factors_tuple_index<FACTOR_CONTAINER_TYPE>();
   std::get<factor_idx>(factors_tuple_).push_back(f);
   return f;
}

template<typename FMC>
inline void factors_storage<FMC>::add_factor_relation(FactorTypeAdapter* f1, FactorTypeAdapter* f2)
{
   add_forward_pass_factor_relation(f1,f2);
   add_backward_pass_factor_relation(f2,f1);
}
template<typename FMC>
void factors_storage<FMC>::add_forward_pass_factor_relation(FactorTypeAdapter* f1, FactorTypeAdapter* f2)
{
   forward_pass_factor_ordering_.clear();
   forward_pass_factor_update_ordering_.clear();
   forward_pass_factor_relation_.push_back({f1,f2});
}
template<typename FMC>
void factors_storage<FMC>::add_backward_pass_factor_relation(FactorTypeAdapter* f1, FactorTypeAdapter* f2)
{
   backward_pass_factor_ordering_.clear();
   backward_pass_factor_update_ordering_.clear();
   backward_pass_factor_relation_.push_back({f1,f2});
}

template<typename FMC>
std::tuple<std::vector<FactorTypeAdapter*>&, std::vector<FactorTypeAdapter*>&> 
factors_storage<FMC>::get_sorted_factors(const Direction d)
{
   assert(number_of_factors() > 0);

   auto get_sorted_factors_impl = [&](std::vector<FactorTypeAdapter*>& ordering, std::vector<FactorTypeAdapter*>& update_ordering, const std::vector<std::array<FactorTypeAdapter*,2>>& factor_rel)
   {
      if(ordering.size() < number_of_factors()) {
         std::tie(ordering, update_ordering) = sort_factors(factor_rel); 
      }
      return std::tie(ordering, update_ordering);
   };

   if(d == Direction::forward) {
      return get_sorted_factors_impl(forward_pass_factor_ordering_, forward_pass_factor_update_ordering_, forward_pass_factor_relation_);
   } else {
      assert(d == Direction::backward);
      return get_sorted_factors_impl(backward_pass_factor_ordering_, backward_pass_factor_update_ordering_, backward_pass_factor_relation_);
   }
}

template<typename FMC>
template<typename FUNC>
void factors_storage<FMC>::for_each_factor(FUNC func) const
{
   for_each_tuple(factors_tuple_, [&func](auto& v) {
         for(auto* f : v) func(*f);
    }); 
}

template<typename FMC>
template<typename FUNC>
void factors_storage<FMC>::for_each_factor(FUNC func)
{
   for_each_tuple(factors_tuple_, [&func](auto& v) {
         for(auto* f : v) func(*f);
    }); 
}

template<typename FMC>
std::size_t factors_storage<FMC>::get_factor_index(const FactorTypeAdapter* f) const
{
   assert(factor_address_to_index_.count(f) == 1);
   return factor_address_to_index_.find(f)->second;
}

template<typename FMC>
template<typename FACTOR_CONTAINER_TYPE>
constexpr std::size_t factors_storage<FMC>::factors_tuple_index()
{
   constexpr auto n = meta::find_index<typename FMC::FactorList, FACTOR_CONTAINER_TYPE>::value;
   static_assert(n < meta::size<typename FMC::FactorList>::value);
   return n;
} 

template<typename FMC>
std::tuple<std::vector<FactorTypeAdapter*>, std::vector<FactorTypeAdapter*>>
factors_storage<FMC>::sort_factors(const std::vector<std::array<FactorTypeAdapter*,2>>& factor_rel)
{

  // assume that factorRel_ describe a DAG. Compute topological sorting
  Topological_Sort::Graph g(number_of_factors());

  for(const auto& rel : factor_rel) {
    assert(factor_address_to_index_.find(rel[0]) != factor_address_to_index_.end());
    assert(factor_address_to_index_.find(rel[1]) != factor_address_to_index_.end());
    const auto f1 = factor_address_to_index_[rel[0]];
    const auto f2 = factor_address_to_index_[rel[1]];
    g.addEdge(f1,f2);
  }

  if(debug()) { std::cout << "sort " << number_of_factors() << " factors subject to " << std::distance(factor_rel.begin(), factor_rel.end()) << " ordering constraints\n"; }

  auto f_sorted = g.topologicalSort();
  assert(f_sorted.size() == number_of_factors());

  std::vector<FactorTypeAdapter*> ordering; 
  ordering.reserve(number_of_factors());
  std::vector<FactorTypeAdapter*> update_ordering;
  for(std::size_t i=0; i<f_sorted.size(); i++) {
     auto* f = factors_[f_sorted[i]];
     ordering.push_back(f);
     if(f->FactorUpdated()) {
        update_ordering.push_back(f);
     }
  }

  return {ordering, update_ordering};
}

} // namespace LPMP

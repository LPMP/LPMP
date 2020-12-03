#pragma once

#include "LP.h"
#include "solver.hxx"
#include "cycle_inequalities.hxx"
#include "tree_decomposition.hxx"
#include "arboricity.h"
#include "mrf_input.h"

#include <string>
#include <vector>
#include <functional>

namespace LPMP {

// expects simplex factor as unary and pairwise factors and marg message such that unary factor is on the left side and pairwise factor is on the right side
// possibly use static inheritance instead of virtual functions
template<class FACTOR_MESSAGE_CONNECTION, std::size_t UNARY_FACTOR_NO, std::size_t PAIRWISE_FACTOR_NO, std::size_t LEFT_MESSAGE_NO, std::size_t RIGHT_MESSAGE_NO>
class mrf_constructor {
public:
   using FMC = FACTOR_MESSAGE_CONNECTION;

   using UnaryFactorContainer = meta::at_c<typename FMC::FactorList, UNARY_FACTOR_NO>;
   using UnaryFactorType = typename UnaryFactorContainer::FactorType;
   using PairwiseFactorContainer = meta::at_c<typename FMC::FactorList, PAIRWISE_FACTOR_NO>;
   using PairwiseFactorType = typename PairwiseFactorContainer::FactorType;
   using LeftMessageContainer = typename meta::at_c<typename FMC::MessageList, LEFT_MESSAGE_NO>::MessageContainerType;
   using LeftMessageType = typename LeftMessageContainer::MessageType;
   using RightMessageContainer = typename meta::at_c<typename FMC::MessageList, RIGHT_MESSAGE_NO>::MessageContainerType;
   using RightMessageType = typename RightMessageContainer::MessageType;


   template<typename SOLVER>
   mrf_constructor(SOLVER& solver) : lp_(&solver.GetLP()) {}

   mrf_constructor(LP<FMC>* lp) : lp_(lp) {}

   std::vector<FactorTypeAdapter*> get_factors()
   {
      std::vector<FactorTypeAdapter*> factors;
      factors.reserve(get_number_of_variables() + get_number_of_pairwise_factors());
      for(auto* f : unaryFactor_) { factors.push_back(f); }
      for(auto* f : pairwiseFactor_) { factors.push_back(f); }
      return factors;
   }

   template<typename T>
   UnaryFactorContainer* add_unary_factor(const std::initializer_list<T> c)
   {
       return add_unary_factor(c.begin(), c.end());
   }

   template<typename ITERATOR>
   UnaryFactorContainer* add_unary_factor(ITERATOR cost_begin, ITERATOR cost_end)
   {
      auto* u = lp_->template add_factor<UnaryFactorContainer>(cost_begin, cost_end);
      unaryFactor_.push_back(u);
      return u; 

   }

   template<typename VECTOR>
   UnaryFactorContainer* add_unary_factor(const VECTOR& cost)
   {
       return add_unary_factor(cost.begin(), cost.end());
   }
   
   // unary factor was created elsewhere, let mrf know it
   void RegisterUnaryFactor(const std::size_t node_number, UnaryFactorContainer* u)
   {
       assert(false);
      if(node_number >= unaryFactor_.size()) {
         unaryFactor_.resize(node_number+1,nullptr);
      } else {
         if(unaryFactor_[node_number] != nullptr) { throw std::runtime_error("unary factor " + std::to_string(node_number) + " already present"); }
      }
      unaryFactor_[node_number] = u;
   }

   template<typename MATRIX>
   PairwiseFactorContainer* add_pairwise_factor(const std::size_t var1, const std::size_t var2, const MATRIX& cost)
   { 
      auto* p = add_pairwise_factor(var1, var2);
      for(std::size_t i=0; i<get_number_of_labels(var1); ++i) {
          for(std::size_t j=0; j<get_number_of_labels(var2); ++j) {
              p->get_factor()->cost(i,j) = cost(i,j);
          }
      }
      return p;
   }

   PairwiseFactorContainer* add_pairwise_factor(const std::size_t var1, const std::size_t var2)
   { 
      assert(var1<var2);
      assert(!has_pairwise_factor(var1,var2));
      auto* p = lp_->template add_factor<PairwiseFactorContainer>(get_number_of_labels(var1), get_number_of_labels(var2));
      pairwiseFactor_.push_back(p);
      pairwiseIndices_.push_back(std::array<std::size_t,2>({var1,var2}));
      const std::size_t factorId = pairwiseFactor_.size()-1;
      pairwiseMap_.insert(std::make_pair(std::make_tuple(var1,var2), factorId));

      lp_->template add_message<LeftMessageContainer>(this->get_unary_factor(var1), p);
      lp_->template add_message<RightMessageContainer>(this->get_unary_factor(var2), p);

      return p;
   }

   // TODO: remove this function
   PairwiseFactorContainer* add_empty_pairwise_factor(const std::size_t var1, const std::size_t var2)
   {
      assert(this->pairwiseMap_.find(std::make_tuple(var1,var2)) == this->pairwiseMap_.end()); 
      return this->add_pairwise_factor(var1, var2);
   } 

   UnaryFactorContainer* get_unary_factor(const std::size_t i) const { assert(i<unaryFactor_.size()); return unaryFactor_[i]; }
   PairwiseFactorContainer* get_pairwise_factor(const std::size_t i) const { assert(i<pairwiseFactor_.size()); return pairwiseFactor_[i]; }
   PairwiseFactorContainer* get_pairwise_factor(const std::size_t i, const std::size_t j) const { 
      assert(i<j);    
      assert(j<unaryFactor_.size());
      const std::size_t factor_id = get_pairwise_factor_id(i,j);
      return pairwiseFactor_[factor_id]; 
   }

   std::pair<LeftMessageContainer*, RightMessageContainer*> get_unary_messages(PairwiseFactorContainer* p) const
   {
      return {get_left_message(p), get_right_message(p)};
   }

   std::pair<LeftMessageContainer*,RightMessageContainer*> get_unary_messages(const std::size_t i, const std::size_t j) const
   {
      assert(i < j);
      auto* p = get_pairwise_factor(i,j);
      return get_unary_messages(p);
   }

   LeftMessageContainer* get_left_message(PairwiseFactorContainer* p) const
   {
      auto msgs = p->template get_messages<LeftMessageContainer>();
      assert(msgs.size() == 1);
      return msgs[0];
   }

   RightMessageContainer* get_right_message(PairwiseFactorContainer* p) const
   {
      auto msgs = p->template get_messages<RightMessageContainer>();
      assert(msgs.size() == 1);
      return msgs[0];
   }

   std::size_t get_number_of_variables() const 
   { 
      assert(!(unaryFactor_.size() == 0 && pairwiseFactor_.size() > 0)); 
      return unaryFactor_.size(); 
   } // note: this is not a good idea, if unaryFactors are populated elsewhere: take maximum in pairwise factor indices then.
   std::size_t get_pairwise_factor_id(const std::size_t var1, const std::size_t var2) const 
   {
      assert(var1<var2);
      assert(pairwiseMap_.find(std::make_tuple(var1,var2)) != pairwiseMap_.end());
      return pairwiseMap_.find(std::make_tuple(var1,var2))->second; 
   }
   bool has_pairwise_factor(const std::size_t var1, const std::size_t var2) const
   {
      assert(var1<var2);
      if(pairwiseMap_.find(std::make_tuple(var1,var2)) != pairwiseMap_.end()) { 
         return true;
      } else {
         return false;
      }
   }
   std::size_t get_number_of_pairwise_factors() const { return pairwiseFactor_.size(); }
   std::array<std::size_t,2> get_pairwise_variables(const std::size_t factorNo) const { assert(factorNo < pairwiseIndices_.size()); return pairwiseIndices_[factorNo]; }
   std::size_t get_number_of_labels(const std::size_t i) const { assert(i < unaryFactor_.size()); return unaryFactor_[i]->get_factor()->size(); }
   REAL GetPairwiseValue(const std::size_t factorId, const std::size_t i1, const std::size_t i2) const
   {
      assert(i1 < get_number_of_labels( get_pairwise_variables(factorId)[0] ));
      assert(i2 < get_number_of_labels( get_pairwise_variables(factorId)[1] ));
      return (*pairwiseFactor_[factorId]->get_factor())(i1,i2);
   }

   void order_factors()
   {
      if(unaryFactor_.size() > 1) {
         for(auto it=unaryFactor_.begin(); it+1 != unaryFactor_.end(); ++it) {
            assert(*it != nullptr);
            lp_->add_factor_relation(*it,*(it+1));
         }
      }

      for(std::size_t p=0; p<get_number_of_pairwise_factors(); ++p) {
          const auto [i,j] = get_pairwise_variables(p);
          assert(i<j);
          auto* f_i = get_unary_factor(i);
          auto* f_j = get_unary_factor(j);
          auto* f_p = get_pairwise_factor(p);
          lp_->add_factor_relation(f_i, f_p);
          lp_->add_factor_relation(f_p, f_j);
      } 
   }

   template<typename SOLVER>
   void Construct(SOLVER& pd) 
   {
      if(debug()) { std::cout << "Construct MRF problem with " << unaryFactor_.size() << " unary factors and " << pairwiseFactor_.size() << " pairwise factors\n"; }

      order_factors();

      assert(pairwiseIndices_.size() == pairwiseFactor_.size());
   }

   std::size_t Tighten(const std::size_t no_constraints_to_add)
   {
      return 0;
   }

   template<typename STREAM>
   void WritePrimal(STREAM& s) const 
   {
      if(unaryFactor_.size() > 0) {
         for(std::size_t i=0; i<unaryFactor_.size()-1; ++i) {
            s << unaryFactor_[i]->get_factor()->primal() << ", ";
         }
         auto* f = unaryFactor_.back();
         s << f->get_factor()->primal() << "\n";
      }
   }

  // build tree of unary and pairwise factors
  factor_tree<FMC> add_tree(std::vector<PairwiseFactorContainer*> p)
  {
     factor_tree<FMC> t;
     assert(p.size() > 0);

     // extract messages joining unaries and pairwise
     std::set<UnaryFactorContainer*> u_set; // u_set is not strictly needed, but helps in checking correctness of algorithm
     std::set<PairwiseFactorContainer*> p_set;
     UnaryFactorContainer* root = nullptr;
     for(auto* f : p) {
        p_set.insert(f); 

        auto left_msgs = f->template get_messages<LeftMessageContainer>();
        assert(left_msgs.size() == 1);
        UnaryFactorContainer* left = (*left_msgs.begin())->GetLeftFactor();
        u_set.insert(left);

        auto right_msgs = f->template get_messages<RightMessageContainer>();
        assert(right_msgs.size() == 1);
        UnaryFactorContainer* right = (*right_msgs.begin())->GetLeftFactor(); 
        u_set.insert(right);

        assert(left != right);

        root = left;
     }

     // assume some arbitrary root. build tree recursively from root
     assert(u_set.size() == p_set.size() + 1);

     std::deque<UnaryFactorContainer*> u_stack;
     u_stack.push_back(root);
     std::vector<std::tuple<Chirality, PairwiseFactorContainer*>> ordered_msgs;
     ordered_msgs.reserve(2*p.size());

     while(!u_stack.empty()) {
        auto* f = u_stack.front();
        u_stack.pop_front();
        u_set.erase(f);

        {
           auto msgs = f->template get_messages<LeftMessageContainer>();
           for(auto it=msgs.begin(); it!= msgs.end(); ++it) {
              auto* p_cand = (*it)->GetRightFactor();
              if(p_set.find(p_cand) != p_set.end()) {
                 p_set.erase(p_cand);
                 ordered_msgs.push_back(std::make_tuple(Chirality::left, p_cand));
                 // search for the other unary connected to p_cand
                 auto msgs_other = p_cand->template get_messages<RightMessageContainer>();
                 assert(msgs_other.size() == 1);
                 auto* u_other = msgs_other[0]->GetLeftFactor();
                 assert(u_set.find(u_other) != u_set.end());
                 //t.AddMessage(*(msgs_other.begin()), Chirality::right);
                 u_stack.push_back(u_other);
              }
           }
        }

        {
           auto msgs = f->template get_messages<RightMessageContainer>();
           for(auto it=msgs.begin(); it!= msgs.end(); ++it) {
              auto* p_cand = (*it)->GetRightFactor();
              if(p_set.find(p_cand) != p_set.end()) {
                 p_set.erase(p_cand);
                 //t.AddMessage((*it), Chirality::left); // or right?
                 ordered_msgs.push_back(std::make_tuple(Chirality::right, p_cand));
                 // search for the other unary connected to p_cand
                 auto msgs_other = p_cand->template get_messages<LeftMessageContainer>();
                 assert(msgs_other.size() == 1);
                 auto* u_other = msgs_other[0]->GetLeftFactor();
                 assert(u_set.find(u_other) != u_set.end());
                 //t.AddMessage(*(msgs_other.begin()), Chirality::right);
                 u_stack.push_back(u_other);
              }
           }
        }
     }

     assert(p_set.empty());

     std::reverse(ordered_msgs.begin(), ordered_msgs.end());
     for(auto e : ordered_msgs) {
         auto* p = std::get<1>(e);
         auto left_msgs = p->template get_messages<LeftMessageContainer>();
         assert(left_msgs.size() == 1);
         auto right_msgs = p->template get_messages<RightMessageContainer>();
         assert(right_msgs.size() == 1);

         if(std::get<0>(e) == Chirality::left) {
             t.add_message(right_msgs[0], Chirality::right);
             t.add_message(left_msgs[0], Chirality::left);
         } else {
             assert(std::get<0>(e) == Chirality::right);
             t.add_message(left_msgs[0], Chirality::right);
             t.add_message(right_msgs[0], Chirality::left);
         }
     }

     return t;
  }

  // compute forest cover of MRF and add each resulting tree

  auto compute_forest_cover() { return compute_forest_cover(pairwiseIndices_); }

  std::vector<factor_tree<FMC>> compute_forest_cover(const std::vector<std::array<std::size_t,2>>& pairwiseIndices)
  {
     UndirectedGraph g(unaryFactor_.size(), pairwiseFactor_.size());
     for(auto e : pairwiseIndices) {
        const auto i = e[0];
        const auto j = e[1];
        g.AddEdge(i,j,1);
     }
     
     const std::size_t forest_num = g.Solve();
     if(diagnostics()) { std::cout << "decomposed mrf into " << forest_num << " trees\n"; }
     std::vector<factor_tree<FMC>> trees;

     std::vector<int> parents(unaryFactor_.size());
     for(std::size_t k=0; k<forest_num; ++k) {
        g.GetForestParents(k, parents.data()); // possible GetForestEdges will return pairwise ids, hence get_pairwise_factor(id) then can be called, which is faster

        union_find uf(unaryFactor_.size());
        for(std::size_t i=0; i<parents.size(); ++i) {
           if(parents[i] != -1) {
              uf.merge(i, parents[i]);
           }
        }
        auto contiguous_ids = uf.get_contiguous_ids();
        
        const std::size_t no_trees = uf.count();
        std::vector<std::vector<PairwiseFactorContainer*>> pairwise(no_trees);

        for(std::size_t i=0; i<parents.size(); ++i) {
           if(parents[i] != -1) {
              const std::size_t tree_id = contiguous_ids[uf.find(i)];
              const std::size_t j = parents[i];
              assert(tree_id == contiguous_ids[uf.find(j)]);
              PairwiseFactorContainer* p = get_pairwise_factor(std::min(i,j), std::max(i,j));
              pairwise[tree_id].push_back(p);
           } 
        }
        for(std::size_t t=0; t<no_trees; ++t) {
           if(pairwise[t].size() > 0) {
              trees.push_back(add_tree(pairwise[t])); 
           }
        }
     }

     auto check_pairwise_factors_present = [&trees]() -> std::size_t {
        std::size_t pairwise_factors_in_trees = 0;
        for(auto& tree : trees) {
           for(auto* f : tree.factors_) {
              if(dynamic_cast<PairwiseFactorContainer*>(f)) {
                 pairwise_factors_in_trees++;
              }
           }
        }
        return pairwise_factors_in_trees;
     };

     //assert(check_pairwise_factors_present() == pairwiseIndices.size());

     return std::move(trees);
  }

  void send_messages_to_unaries()
  {
     for (auto *p : pairwiseFactor_)
     {
        auto msgs = get_unary_messages(p);
        std::get<0>(msgs)->send_message_to_right(1.0);
        std::get<1>(msgs)->send_message_to_right(1.0);

        std::get<0>(msgs)->send_message_to_left(0.5);
        std::get<1>(msgs)->send_message_to_left(1.0);
        std::get<0>(msgs)->send_message_to_left(1.0);
     }
  }

  void construct(const mrf_input& input)
  {
      // first input the unaries, as pairwise potentials need them to be able to link to them
      for(std::size_t i=0; i<input.no_variables(); ++i) {
          this->add_unary_factor(input.unaries[i]);
      }

      assert(input.pairwise_indices.size() == input.pairwise_values.dim1());
      for(std::size_t i=0; i<input.pairwise_indices.size(); ++i) {
          const auto var1 = input.pairwise_indices[i][0];
          const auto var2 = input.pairwise_indices[i][1];
          auto pairwise_cost = input.pairwise_values[i];
          this->add_pairwise_factor(var1,var2, pairwise_cost);
      }
  }

  LP<FMC>* get_lp() const { return lp_; }

protected:
   std::vector<UnaryFactorContainer*> unaryFactor_;
   std::vector<PairwiseFactorContainer*> pairwiseFactor_;
   
   std::vector<std::array<std::size_t,2>> pairwiseIndices_;

   std::map<std::tuple<std::size_t,std::size_t>, std::size_t> pairwiseMap_; // given two sorted indices, return factorId belonging to that index.

   //std::size_t unaryFactorIndexBegin_, unaryFactorIndexEnd_; // do zrobienia: not needed anymore

   LP<FMC>* lp_;
};

// derives from a given mrf problem constructor and adds tightening capabilities on top of it, as implemented in cycle_inequalities and proposed by David Sontag
template<class MRF_PROBLEM_CONSTRUCTOR,
   std::size_t TERNARY_FACTOR_NO, std::size_t PAIRWISE_TRIPLET_MESSAGE12_NO, std::size_t PAIRWISE_TRIPLET_MESSAGE13_NO, std::size_t PAIRWISE_TRIPLET_MESSAGE23_NO> // the last indices indicate triplet factor and possible messages
class tightening_mrf_constructor : public MRF_PROBLEM_CONSTRUCTOR
{
public:
   using MRFPC = MRF_PROBLEM_CONSTRUCTOR;
   using FMC = typename MRFPC::FMC;
   using MrfConstructorType = tightening_mrf_constructor<MRFPC, TERNARY_FACTOR_NO, PAIRWISE_TRIPLET_MESSAGE12_NO, PAIRWISE_TRIPLET_MESSAGE13_NO, PAIRWISE_TRIPLET_MESSAGE23_NO>;

   using TripletFactorContainer = meta::at_c<typename FMC::FactorList, TERNARY_FACTOR_NO>; 
   using TripletFactor = typename TripletFactorContainer::FactorType;
   using PairwiseTripletMessage12Container = typename meta::at_c<typename FMC::MessageList, PAIRWISE_TRIPLET_MESSAGE12_NO>::MessageContainerType;
   using PairwiseTripletMessage13Container = typename meta::at_c<typename FMC::MessageList, PAIRWISE_TRIPLET_MESSAGE13_NO>::MessageContainerType;
   using PairwiseTripletMessage23Container = typename meta::at_c<typename FMC::MessageList, PAIRWISE_TRIPLET_MESSAGE23_NO>::MessageContainerType;

   template<typename SOLVER>
   tightening_mrf_constructor(SOLVER& pd)
      : MRF_PROBLEM_CONSTRUCTOR(pd)
   {}

   TripletFactorContainer* add_triplet_factor(const std::size_t var1, const std::size_t var2, const std::size_t var3, const std::vector<REAL>& cost)
   {
      assert(var1<var2 && var2<var3);
      assert(var3<this->get_number_of_variables());
      assert(tripletMap_.find(std::array<std::size_t,3>({var1,var2,var3})) == tripletMap_.end());
      
      assert(this->pairwiseMap_.find(std::make_tuple(var1,var2)) != this->pairwiseMap_.end());
      assert(this->pairwiseMap_.find(std::make_tuple(var1,var3)) != this->pairwiseMap_.end());
      assert(this->pairwiseMap_.find(std::make_tuple(var2,var3)) != this->pairwiseMap_.end());

      const std::size_t factor12Id = this->pairwiseMap_.find(std::make_tuple(var1,var2))->second;
      const std::size_t factor13Id = this->pairwiseMap_.find(std::make_tuple(var1,var3))->second;
      const std::size_t factor23Id = this->pairwiseMap_.find(std::make_tuple(var2,var3))->second;

      TripletFactorContainer* t = this->lp_->template add_factor<TripletFactorContainer>(this->get_number_of_labels(var1), this->get_number_of_labels(var2), this->get_number_of_labels(var3));
      tripletFactor_.push_back(t);
      tripletIndices_.push_back(std::array<std::size_t,3>({var1,var2,var3}));
      const std::size_t factorId = tripletFactor_.size()-1;
      tripletMap_.insert(std::make_pair(std::array<std::size_t,3>({var1,var2,var3}), factorId));

      LinkPairwiseTripletFactor<PairwiseTripletMessage12Container>(factor12Id,factorId);
      LinkPairwiseTripletFactor<PairwiseTripletMessage13Container>(factor13Id,factorId);
      LinkPairwiseTripletFactor<PairwiseTripletMessage23Container>(factor23Id,factorId); 

      //this->lp_->ForwardPassFactorRelation(this->get_pairwise_factor(factor12Id),t);
      //this->lp_->ForwardPassFactorRelation(this->get_pairwise_factor(factor13Id),t);
      //this->lp_->ForwardPassFactorRelation(t,this->get_pairwise_factor(factor23Id));

      //this->lp_->BackwardPassFactorRelation(this->get_pairwise_factor(factor23Id),t);
      //this->lp_->BackwardPassFactorRelation(this->get_pairwise_factor(factor13Id),t);
      //this->lp_->BackwardPassFactorRelation(t,this->get_pairwise_factor(factor12Id));
      return t;
   }
   template<typename PAIRWISE_TRIPLET_MESSAGE_CONTAINER>
   void LinkPairwiseTripletFactor(const std::size_t pairwiseFactorId, const std::size_t tripletFactorId)
   {
      using PairwiseTripletMessageType = typename PAIRWISE_TRIPLET_MESSAGE_CONTAINER::MessageType;

      typename MRFPC::PairwiseFactorContainer* const p = this->pairwiseFactor_[pairwiseFactorId];
      assert(this->pairwiseIndices_[pairwiseFactorId][0] < this->pairwiseIndices_[pairwiseFactorId][1]);

      TripletFactorContainer* const t = tripletFactor_[tripletFactorId];
      const std::size_t tripletVar1 = tripletIndices_[tripletFactorId][0];
      const std::size_t tripletVar2 = tripletIndices_[tripletFactorId][1];
      const std::size_t tripletVar3 = tripletIndices_[tripletFactorId][2];
      assert(tripletVar1 < tripletVar2 && tripletVar2 < tripletVar3);
      const std::size_t tripletDim1 = this->get_number_of_labels(tripletVar1);
      const std::size_t tripletDim2 = this->get_number_of_labels(tripletVar2);
      const std::size_t tripletDim3 = this->get_number_of_labels(tripletVar3);
         
      assert(this->get_number_of_labels(this->pairwiseIndices_[pairwiseFactorId][0]) * this->get_number_of_labels(this->pairwiseIndices_[pairwiseFactorId][1]) == p->get_factor()->size());

      this->lp_->template add_message<PAIRWISE_TRIPLET_MESSAGE_CONTAINER>(p, t, tripletDim1, tripletDim2, tripletDim3);
   }
   std::size_t get_number_of_triplet_factors() const { return tripletFactor_.size(); }

   std::array<std::size_t,3> get_triplet_indices(const std::size_t factor_id)
   {
      assert(factor_id < get_number_of_triplet_factors());
      return tripletIndices_[factor_id];
   }

   // do zrobienia: use references for pi
   bool add_tightening_triplet(const std::size_t var1, const std::size_t var2, const std::size_t var3)//, const std::vector<std::size_t> pi1, const std::vector<std::size_t> pi2, const std::vector<std::size_t> pi3)
   {
      assert(var1 < var2 && var2 < var3 && var3 < this->get_number_of_variables());
      if(tripletMap_.count(std::array<std::size_t,3>({var1,var2,var3})) == 0) {
         // first check whether necessary pairwise factors are present. If not, add them.
         if(this->pairwiseMap_.find(std::make_tuple(var1,var2)) == this->pairwiseMap_.end()) {
            this->add_empty_pairwise_factor(var1,var2);
         }
         if(this->pairwiseMap_.find(std::make_tuple(var1,var3)) == this->pairwiseMap_.end()) {
            this->add_empty_pairwise_factor(var1,var3);
         }
         if(this->pairwiseMap_.find(std::make_tuple(var2,var3)) == this->pairwiseMap_.end()) {
            this->add_empty_pairwise_factor(var2,var3);
         }

         const auto tripletSize = this->get_number_of_labels(var1) *  this->get_number_of_labels(var2) * this->get_number_of_labels(var3);
         add_triplet_factor(var1,var2,var3, std::vector<REAL>(tripletSize,0.0));
         return true;
      } else {
         return false;
      }
   }

   std::size_t add_triplets(const std::vector<triplet_candidate>& tc, const std::size_t max_triplets_to_add = std::numeric_limits<std::size_t>::max())
   {
      std::size_t no_triplets_added = 0;
      for(const auto t : tc) {
         if(add_tightening_triplet(t.i, t.j, t.k)) {
            ++no_triplets_added;
            if(no_triplets_added >= max_triplets_to_add) {
               break;
            }
         }
      }
      return no_triplets_added; 
   }

   std::size_t Tighten(const std::size_t noTripletsToAdd)
   {
      assert(noTripletsToAdd > 0);
      if(debug()) {
         std::cout << "Tighten mrf with cycle inequalities, no triplets to add = " << noTripletsToAdd << "\n";
      }

      //auto fp = [this](const std::size_t v1, const std::size_t v2, const std::size_t v3) { return this->add_tightening_triplet(v1,v2,v3); }; // do zrobienia: do not give this via template, as Cycle already has gm_ object.

      triplet_search<typename std::remove_reference<decltype(*this)>::type> triplets(*this, eps);
      if(debug()) { std::cout << "search for triplets\n"; }
      auto triplet_candidates = triplets.search();
      if(debug()) { std::cout << "done\n"; }
      std::size_t no_triplets_added = add_triplets(triplet_candidates, noTripletsToAdd);
      if(diagnostics()) { std::cout << "added " << no_triplets_added << " by triplet search\n"; }

      if(no_triplets_added < 0.2*noTripletsToAdd) {
         if(debug()) {
            std::cout << "----------------- do cycle search with k-projection graph---------------\n";
         }
         k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(*this)>::type, false> cycle_search(*this, eps);
         triplet_candidates = cycle_search.search();
         if(debug()) { std::cout << "... done\n"; }
         const std::size_t no_triplet_k_projection_graph = add_triplets(triplet_candidates, noTripletsToAdd-no_triplets_added);
         if(diagnostics()) {std::cout << "added " << no_triplet_k_projection_graph << " by cycle search in k-projection graph\n"; }
         no_triplets_added += no_triplet_k_projection_graph;

         // search in expanded projection graph
         if(no_triplets_added < 0.2*noTripletsToAdd) {
            if(debug()) {
               std::cout << "----------------- do cycle search with expanded projection graph---------------\n";
            }
            k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(*this)>::type, true> cycle_search(*this, eps);
            triplet_candidates = cycle_search.search();
            if(debug()) { std::cout << "... done\n"; }
            const std::size_t no_triplet_k_projection_graph = add_triplets(triplet_candidates, noTripletsToAdd-no_triplets_added);
            if(diagnostics()) { std::cout << "added " << no_triplet_k_projection_graph << " by cycle search in expanded projection graph\n";
            }
            no_triplets_added += no_triplet_k_projection_graph;
         }
      }

      if(debug()) {
         std::cout << "no triplets found, now search for triplets with guaranteed increase also smaller than " << eps << "\n";
      }
       
      assert(eps >= std::numeric_limits<REAL>::epsilon());

      if(no_triplets_added == 0) {
         triplet_search<typename std::remove_reference<decltype(*this)>::type> triplets(*this, std::numeric_limits<REAL>::epsilon());
         if(debug()) {
            std::cout << "search for triplets (any will do)\n";
         }
         auto triplet_candidates = triplets.search();
         if(debug()) {
            std::cout << "done\n";
         }
         no_triplets_added = add_triplets(triplet_candidates, noTripletsToAdd);
         if(diagnostics()) {
            std::cout << "added " << no_triplets_added << " by triplet search with no increase\n";
         }
      }

      if(no_triplets_added == 0) {
         if(debug()) {
            std::cout << "----------------- do cycle search with k-projection graph, add all---------------\n";
         }
         k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(*this)>::type, false> cycle_search(*this, std::numeric_limits<REAL>::epsilon());
         triplet_candidates = cycle_search.search();
         if(debug()) { std::cout << "... done\n"; }
         const std::size_t no_triplet_k_projection_graph = add_triplets(triplet_candidates, noTripletsToAdd-no_triplets_added);
         if(diagnostics()) {
            std::cout << "added " << no_triplet_k_projection_graph << " by cycle search in k-projection grapa with no increase\n";
         }
         no_triplets_added += no_triplet_k_projection_graph;
      }

      if(no_triplets_added == 0) {
         if(debug()) {
            std::cout << "----------------- do cycle search with expanded projection graph, add all---------------\n";
         }
         k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(*this)>::type, true> cycle_search(*this, std::numeric_limits<REAL>::epsilon());
         triplet_candidates = cycle_search.search();
         if(debug()) { std::cout << "... done\n"; }
         const std::size_t no_triplet_k_projection_graph = add_triplets(triplet_candidates, noTripletsToAdd-no_triplets_added);
         if(diagnostics()) {
            std::cout << "added " << no_triplet_k_projection_graph << " by cycle search in expanded projection graph with no increase\n";
         }
         no_triplets_added += no_triplet_k_projection_graph;
      }

      return no_triplets_added;
   }

  void send_messages_to_unaries()
  {
     for(auto* t : tripletFactor_) {
        auto msgs = get_triplet_to_pairwise_messages(t);
        std::get<0>(msgs)->send_message_to_right(1.0);
        std::get<1>(msgs)->send_message_to_right(1.0);
        std::get<2>(msgs)->send_message_to_right(1.0);

        std::get<0>(msgs)->send_message_to_left(1.0/3.0);
        std::get<1>(msgs)->send_message_to_left(1.0/3.0);
        std::get<2>(msgs)->send_message_to_left(1.0);
        std::get<1>(msgs)->send_message_to_left(0.5);
        std::get<0>(msgs)->send_message_to_left(1.0); 
     }

     MRF_PROBLEM_CONSTRUCTOR::send_messages_to_unaries();
  } 

protected:

  std::tuple< PairwiseTripletMessage12Container*, PairwiseTripletMessage13Container*, PairwiseTripletMessage23Container* >
  get_triplet_to_pairwise_messages(TripletFactorContainer* t)
  {
     auto msgs12 = t->template get_messages<PairwiseTripletMessage12Container>();
     auto msgs13 = t->template get_messages<PairwiseTripletMessage13Container>();
     auto msgs23 = t->template get_messages<PairwiseTripletMessage23Container>();
     assert(msgs12.size() == 1 && msgs13.size() == 1 && msgs23.size() == 1);
     return {msgs12[0], msgs13[0], msgs23[0]};
  }

   std::vector<TripletFactorContainer*> tripletFactor_;
   std::vector<std::array<std::size_t,3>> tripletIndices_;
   std::map<std::array<std::size_t,3>, std::size_t> tripletMap_; // given two sorted indices, return factorId belonging to that index.
};

} // end namespace LPMP

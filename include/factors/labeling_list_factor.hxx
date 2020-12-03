#ifndef LPMP_LABELING_LIST_FACTOR_HXX
#define LPMP_LABELING_LIST_FACTOR_HXX

#include <array>
#include <bitset>
#include <cassert>
#include <iostream>
#include "vector.hxx"
#include "serialization.hxx"

// to do: make _impl functions private
// make functions static whenever they can be

namespace LPMP {

template<std::size_t... LABELS>
struct labeling {

   template<std::size_t LABEL_NO, std::size_t LABEL, std::size_t... LABELS_REST>
   constexpr static
   typename std::enable_if<LABEL_NO == 0,std::size_t>::type label_impl()
   {
     static_assert(LABEL == 0 || LABEL == 1,"");
      return LABEL;
   }

   template<std::size_t LABEL_NO, std::size_t LABEL, std::size_t... LABELS_REST>
   constexpr static
   typename std::enable_if<(LABEL_NO > 0),std::size_t>::type label_impl()
   {
      return label_impl<LABEL_NO-1, LABELS_REST...>();
   }

   template<std::size_t LABEL_NO>
   constexpr static std::size_t label()
   { 
      static_assert(LABEL_NO < sizeof...(LABELS), "label number must be smaller than number of labels");
      return label_impl<LABEL_NO, LABELS...>();
   } 

   constexpr static std::size_t no_labels()
   {
      return sizeof...(LABELS);
   }

   template<std::size_t I, std::size_t... LABELS_REST>
   static typename std::enable_if<(I == sizeof...(LABELS)),bool>::type
   matches_impl(const std::bitset< sizeof...(LABELS)>& l)
   {
      return true;
   }

   template<std::size_t I, std::size_t LABEL, std::size_t... LABELS_REST>
   static typename std::enable_if<(I < sizeof...(LABELS)),bool>::type
   matches_impl(const std::bitset< sizeof...(LABELS)>& l)
   {
      if(l[I] != LABEL) {
         return false;
      } else {
         return matches_impl<I+1, LABELS_REST...>(l);
      }
   }

   static bool matches(const std::bitset< sizeof...(LABELS)>& l)
   {
      return matches_impl<0, LABELS...>(l); 
   }

   template<std::size_t I, std::size_t... LABELS_REST>
   static typename std::enable_if<(I == sizeof...(LABELS))>::type
   get_labeling_impl(const std::bitset< sizeof...(LABELS)>& l)
   {
      return;
   }

   template<std::size_t I, std::size_t LABEL, std::size_t... LABELS_REST>
   constexpr static typename std::enable_if<(I < sizeof...(LABELS))>::type
   get_labeling_impl(std::bitset< sizeof...(LABELS)>& l)
   {
     l[I] = (LABEL == 0) ? false : true;
     return get_labeling_impl<I+1, LABELS_REST...>(l);
   }

   constexpr static std::bitset<no_labels()> get_labeling()
   {
     std::bitset<no_labels()> l;
     get_labeling_impl<0, LABELS...>(l);
     return l;
   }
};

template<typename... LABELINGS> // all labels must be instances of labeling
struct labelings
{
   /* destructor makes labelings a non-constant type
   ~labelings()
   {
      static_assert(sizeof...(LABELINGS) > 0, "at least one labeling must be present");
      // to do: check whether each label occurs at most once and each labeling has same number of labels.
   }
   */

   template<std::size_t LABELING_NO, std::size_t LABEL_NO, typename LABELING, typename... LABELINGS_REST>
   constexpr static typename std::enable_if<LABELING_NO == 0,std::size_t>::type 
   get_label()
   {
      return LABELING::template label<LABEL_NO>();
   }

   template<std::size_t LABELING_NO, std::size_t LABEL_NO, typename LABELING, typename... LABELINGS_REST>
   constexpr static typename std::enable_if<(LABELING_NO > 0),std::size_t>::type 
   get_label()
   {
      return get_label<LABELING_NO-1, LABEL_NO, LABELINGS_REST...>();
   }

   template<std::size_t LABELING_NO, std::size_t LABEL_NO>
   constexpr static std::size_t label()
   {
      static_assert(sizeof...(LABELINGS) > 0, "at least one labeling must be present");
      // to do: check whether each label occurs at most once and each labeling has same number of labels.

      static_assert(LABELING_NO < sizeof...(LABELINGS), "labeling number must be smaller than number of labelings");
      return get_label<LABELING_NO,LABEL_NO,LABELINGS...>();
   }

   template<std::size_t IDX, typename LABELING, typename... LABELINGS_REST>
   constexpr static std::bitset<LABELING::no_labels()> get_labeling_impl()
   {
      if constexpr(IDX == 0)
         return LABELING::get_labeling();
      else
         return get_labeling_impl<IDX-1, LABELINGS_REST...>(); 
   }
   template<std::size_t LABELING_NO>
   constexpr static auto get_labeling()
   {
      return get_labeling_impl<LABELING_NO, LABELINGS...>();
   }

   template<std::size_t LABELING_NO>
   static auto get_labeling_impl(const std::size_t labeling_no)
   {
       if(labeling_no == LABELING_NO)
           return get_labeling<LABELING_NO>();
       if constexpr(LABELING_NO+1 < no_labelings())
           return get_labeling_impl<LABELING_NO+1>(labeling_no);
       else 
           throw std::runtime_error("labeling number larger than number of labels"); // TODO: terminate
   }
   static auto get_labeling(const std::size_t labeling_no)
   {
       assert(labeling_no < no_labelings());
       return get_labeling_impl<0>(labeling_no); 
   }

   constexpr static std::size_t no_labelings()
   {
      return sizeof...(LABELINGS);
   }

   template<typename LABELING, typename... LABELINGS_REST>
   constexpr static std::size_t no_labels_impl()
   {
      return LABELING::no_labels();
   }
   constexpr static std::size_t no_labels()
   {
      return no_labels_impl<LABELINGS...>();
   }

   template<std::size_t I, typename... LABELINGS_REST>
   static typename std::enable_if<(I >= sizeof...(LABELINGS)), std::size_t>::type
   matching_labeling_impl(const std::bitset< no_labels()>& l)
   {
      return no_labelings();
   }

   template<std::size_t I, typename LABELING, typename... LABELINGS_REST>
   static typename std::enable_if<(I < sizeof...(LABELINGS)), std::size_t>::type
   matching_labeling_impl(const std::bitset< no_labels()>& l)
   {
      if(LABELING::matches(l)) {
         return I;
      } else {
         return matching_labeling_impl<I+1,LABELINGS_REST...>(l);
      }
   }

   static std::size_t matching_labeling(const std::bitset< no_labels()>& l)
   {
      return matching_labeling_impl<0, LABELINGS...>(l); 
   }

   template<std::size_t I, typename... LABELINGS_REST>
   static typename std::enable_if<(I >= sizeof...(LABELINGS)), std::bitset<no_labels()>>::type
   labeling_impl(const std::size_t no)
   {
      assert(false);
      return std::bitset<no_labels()>();
   }

   template<std::size_t I, typename LABELING, typename... LABELINGS_REST>
   static typename std::enable_if<(I < sizeof...(LABELINGS)), std::bitset<no_labels()>>::type
   labeling_impl(const std::size_t no)
   {
     if(I == no) {
       return LABELING::get_labeling();
     } else {
       return labeling_impl<I+1,LABELINGS_REST...>(no);
     }
   }

   static std::bitset<no_labels()> labeling(const std::size_t no)
   {
     return labeling_impl<0, LABELINGS...>(no);
   }

};

template<>
class labelings<>
{
   template<std::size_t LABELING_NO, std::size_t LABEL_NO>
   constexpr static std::size_t label()
   {
      return 42;
   }

   constexpr static std::size_t no_labelings()
   {
      return 0;
   }

   constexpr static std::size_t no_labels()
   {
      return 0;
   }

   template<typename VEC>
   static std::size_t matching_labeling(const VEC& l)
   {
      return 0;
   }
};

template<typename LABELINGS, bool IMPLICIT_ORIGIN>
class labeling_factor : public array<double, LABELINGS::no_labelings()>
{
public:
   labeling_factor() 
   {
      std::fill(this->begin(), this->end(), 0.0);
   }
   ~labeling_factor()
   {}

   constexpr static bool has_implicit_origin() { return IMPLICIT_ORIGIN; } // means zero label has cost 0 and is not recorded.

   constexpr static std::size_t size() 
   {
      return LABELINGS::no_labelings();
   }
   
   constexpr static std::size_t primal_size() 
   {
      return LABELINGS::no_labels(); 
   }

   constexpr static LABELINGS labelings() { return LABELINGS{}; }

   double LowerBound() const
   {
      assert(this->min() == *std::min_element(this->begin(), this->end()));
      if(has_implicit_origin()) {
         return std::min(0.0, this->min());
      } else {
         return this->min();
      }
   }

   double EvaluatePrimal() const
   {
      const std::size_t labeling_no = LABELINGS::matching_labeling(primal_);
      if(labeling_no < size()) {
         return (*this)[labeling_no];
      }
      // check for zero labeling
      if(has_implicit_origin() && primal_.count() == 0) {
         return 0.0;
      }
      return std::numeric_limits<double>::infinity();
   }

   void MaximizePotentialAndComputePrimal()
   {
       if(primal_.count() == primal_.size()) {
           primal_.reset();
           const std::size_t opt_sol = std::min_element(this->begin(), this->end()) - this->begin();
           primal_ = LABELINGS::get_labeling(opt_sol);
       }
   }

   // return two possible variable states
   // branch on current primal vs. not current primal
   void branch_left()
   {
      // set cost of label associated not with primal to infinity, i.e. current label should always be taken
      const std::size_t labeling_no = LABELINGS::matching_labeling(primal_);
      assert(labeling_no < this->size());
      for(std::size_t i=0; i<this->size(); ++i) {
         (*this)[i] = std::numeric_limits<double>::infinity(); 
      }
      // also the zero labeling must be forbidden. How to do? IMPLICIT_ORIGIN must be dropped for this.
      assert(false);
      assert(!has_implicit_origin());
   }

   void branch_right()
   {
      // set cost of primal label to infinity
      assert(EvaluatePrimal() < std::numeric_limits<double>::infinity());
      const std::size_t labeling_no = LABELINGS::matching_labeling(primal_);
      assert(labeling_no < this->size());
      (*this)[labeling_no] = std::numeric_limits<double>::infinity();
      assert(LowerBound() < std::numeric_limits<double>::infinity());
   }

   auto& primal() { return primal_; }
   const auto& primal() const { return primal_; }

   void init_primal() { primal_.set(); }
   template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( binary_data<double>(&(*this)[0], size()) ); }
   template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }

   auto export_variables() { return std::tie( *static_cast<array<double, LABELINGS::no_labelings()>*>(this) ); }

   template<typename EXTERNAL_SOLVER, typename VECTOR>
   void construct_constraints(EXTERNAL_SOLVER& s, VECTOR vars) const
   {
      if(has_implicit_origin()) {
	      s.add_at_most_one_constraint(vars.begin(), vars.end());
      } else {
	      s.add_simplex_constraint(vars.begin(), vars.end());
      }
   }

   template<typename EXTERNAL_SOLVER, typename VECTOR>
   void convert_primal(EXTERNAL_SOLVER& s, VECTOR vars)
   {
      primal_.reset();
      primal_ = s.first_active(vars.begin(), vars.end());
   }

   template<std::size_t IDX, typename FUNC>
   void for_each_labeling_impl(FUNC f) const
   {
      if constexpr(IDX >= LABELINGS::no_labelings()) {
         return;
      } else {
         f(LABELINGS::template get_labeling<IDX>(), (*this)[IDX]);
         for_each_labeling_impl<IDX+1>(f);
      }
   }

   template<typename FUNC>
   void for_each_labeling(FUNC f) const
   {
      for_each_labeling_impl<0>(f);
   }

   template<std::size_t IDX, typename FUNC>
   void for_each_labeling_impl(FUNC f)
   {
      if constexpr(IDX >= LABELINGS::no_labelings()) {
         return;
      } else {
         f(LABELINGS::template get_labeling<IDX>(), (*this)[IDX]);
         for_each_labeling_impl<IDX+1>(f);
      }
   }

   template<typename FUNC>
   void for_each_labeling(FUNC f)
   {
      for_each_labeling_impl<0>(f);
   }

private:
   std::bitset<primal_size()> primal_;
};

// we assume that LEFT_LABELING contains sublabelings of RIGHT_LABELING, where INDICES indicate i-th entry of LEFT_LABELING is mapped to INDICES[i]-th entry of RIGHT_LABELINGS
template<typename LEFT_LABELINGS, typename RIGHT_LABELINGS, std::size_t... INDICES>
class labeling_message {
using type = labeling_message<LEFT_LABELINGS, RIGHT_LABELINGS, INDICES...>;

public:
   using msg_val_type = array<double, LEFT_LABELINGS::no_labelings()>;
   static constexpr std::array<std::size_t, sizeof...(INDICES)> indices = {INDICES...};

   ~labeling_message() 
   {
      static_assert(sizeof...(INDICES) == LEFT_LABELINGS::no_labels(), "each left label must be matched to a right one");
   }

   template<typename LEFT_LABELING, typename RIGHT_LABELING, std::size_t LEFT_IDX, std::size_t... I_REST>
   constexpr static typename std::enable_if<(LEFT_IDX >= LEFT_LABELING::no_labels()),bool>::type
   matches_impl()
   {
      return true;
   }
   template<typename LEFT_LABELING, typename RIGHT_LABELING, std::size_t LEFT_IDX, std::size_t I, std::size_t... I_REST>
   constexpr static typename std::enable_if<(LEFT_IDX < LEFT_LABELING::no_labels()),bool>::type
   matches_impl()
   {
      if(LEFT_LABELING::template label<LEFT_IDX>() == RIGHT_LABELING::template label<I>()) {
         return type::matches_impl<LEFT_LABELING, RIGHT_LABELING, LEFT_IDX+1, I_REST...>();
      } else {
         return false;
      }
   }

   template<typename LEFT_LABELING, typename RIGHT_LABELING>
   constexpr static bool matches()
   {
      return matches_impl<LEFT_LABELING, RIGHT_LABELING, 0, INDICES...>();
   }

   template<typename RIGHT_LABELING, std::size_t I, typename... LEFT_LABELINGS_REST>
   constexpr static typename std::enable_if<(I >= LEFT_LABELINGS::no_labelings()),std::size_t >::type
   matching_left_labeling_impl(labelings<LEFT_LABELINGS_REST...>)
   {
      return I; // return 1 + number of left labelings;
   }
   template<typename RIGHT_LABELING, std::size_t I, typename LEFT_LABELING, typename... LEFT_LABELINGS_REST>
   constexpr static typename std::enable_if<(I < LEFT_LABELINGS::no_labelings()),std::size_t>::type
   matching_left_labeling_impl(labelings<LEFT_LABELING, LEFT_LABELINGS_REST...>)
   {
      if(matches<LEFT_LABELING, RIGHT_LABELING>()) {
         return I;
      } else {
         return matching_left_labeling_impl<RIGHT_LABELING, I+1>(labelings<LEFT_LABELINGS_REST...>{});
      }
   }

   // we assume that there is as most one matching left labeling
   template<typename RIGHT_LABELING>
   constexpr static std::size_t matching_left_labeling()
   {
      return matching_left_labeling_impl<RIGHT_LABELING,0>(LEFT_LABELINGS{});
   } 

  template<typename RIGHT_FACTOR, std::size_t I, typename... RIGHT_LABELINGS_REST>
  typename std::enable_if<(I >= RIGHT_LABELINGS::no_labelings())>::type 
  compute_msg_impl(msg_val_type& msg_val, const RIGHT_FACTOR& r, double& min_of_labels_not_taken, labelings<RIGHT_LABELINGS_REST...>)
  {
     return;
  }

  template<typename RIGHT_FACTOR, std::size_t I, typename RIGHT_LABELING, typename... RIGHT_LABELINGS_REST>
  typename std::enable_if<(I < RIGHT_LABELINGS::no_labelings())>::type 
  compute_msg_impl(msg_val_type& msg_val, const RIGHT_FACTOR& r, double& min_of_labels_not_taken, labelings<RIGHT_LABELING, RIGHT_LABELINGS_REST...>)
  {
     constexpr std::size_t left_label_number = matching_left_labeling<RIGHT_LABELING>();
     if(left_label_number < msg_val.size()) {
        msg_val[left_label_number] = std::min(msg_val[left_label_number], r[I]);
     } else {
        assert(left_label_number == msg_val.size());
        min_of_labels_not_taken = std::min(min_of_labels_not_taken, r[I]); 
     }
     compute_msg_impl<RIGHT_FACTOR, I+1, RIGHT_LABELINGS_REST...>(msg_val, r, min_of_labels_not_taken, labelings<RIGHT_LABELINGS_REST...>{});
  }

  template<typename RIGHT_FACTOR>
  void compute_msg(msg_val_type& msg_val, const RIGHT_FACTOR& r)
  {
     double min_of_labels_not_taken;
     if(r.has_implicit_origin()) {
        min_of_labels_not_taken = 0.0;
     } else {
        min_of_labels_not_taken = std::numeric_limits<double>::infinity();
     }
     std::fill(msg_val.begin(), msg_val.end(), std::numeric_limits<double>::infinity());
     compute_msg_impl<RIGHT_FACTOR, 0>(msg_val, r, min_of_labels_not_taken, RIGHT_LABELINGS{});
     // note: this is possibly wrong, if r.has_implicit_origin() is false
     for(auto& v : msg_val) {
        v -= min_of_labels_not_taken;
     }
     for(std::size_t i=0; i<msg_val.size(); ++i) {
        assert(!std::isnan(msg_val[i]));
     }
  }

   template<typename RIGHT_FACTOR, typename MSG, std::size_t I, typename... RIGHT_LABELINGS_REST>
   typename std::enable_if<(I >= RIGHT_LABELINGS::no_labelings())>::type 
   repam_right_impl(RIGHT_FACTOR& r, const MSG& msg, labelings<RIGHT_LABELINGS_REST...>) const
   {
      return;
   }

   template<typename RIGHT_FACTOR, typename MSG, std::size_t I, typename RIGHT_LABELING, typename... RIGHT_LABELINGS_REST>
   typename std::enable_if<(I < RIGHT_LABELINGS::no_labelings())>::type 
   repam_right_impl(RIGHT_FACTOR& r, const MSG& msg, labelings<RIGHT_LABELING, RIGHT_LABELINGS_REST...>) const
   {
     constexpr std::size_t left_label_number = matching_left_labeling<RIGHT_LABELING>();
     if(left_label_number < msg.size()) {
        r[I] += msg[left_label_number];
     }
     repam_right_impl<RIGHT_FACTOR, MSG, I+1, RIGHT_LABELINGS_REST...>(r, msg, labelings<RIGHT_LABELINGS_REST...>{});
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void RepamRight(RIGHT_FACTOR& r, const MSG& msg) const
   {
      // msg has dimension equal to number of left labelings;
      // go over all right labelings, find corresponding left labeling (if there is any) and if so, add msg value
      for(std::size_t i=0; i<LEFT_LABELINGS::no_labelings(); ++i) {
         assert(!std::isnan(msg[i]));
         assert(std::isfinite(msg[i]));
      } 
      repam_right_impl<RIGHT_FACTOR, MSG, 0>(r, msg, RIGHT_LABELINGS{}); 
      for(std::size_t i=0; i<r.size(); ++i) {
         assert(!std::isnan(r[i]));
         assert(std::isfinite(r[i]));
      }
   }

   //template<typename RIGHT_FACTOR>
   //void RepamRight(RIGHT_FACTOR& r, const double msg_val, const std::size_t idx) const
   //{
   //    repam_right
   //}

   template<typename RIGHT_FACTOR, typename MSG>
   void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const double omega)
   {
       for(std::size_t i=0; i<r.size(); ++i) {
           assert(std::isfinite(r[i]));
           assert(!std::isnan(r[i]));
       }
      // msg has dimension equal to number of left labelings;
      // go over all right labelings. Then find left labeling corresponding to it, and compute minimum
      msg_val_type msg_val;
      compute_msg(msg_val, r);
      for(std::size_t i=0; i<msg_val.size(); ++i) {
          assert(std::isfinite(msg_val[i]));
          assert(!std::isnan(msg_val[i]));
      } 
      msg -= omega*msg_val;
   }

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& l, const MSG& msg) const
   {
      for(std::size_t i=0; i<LEFT_LABELINGS::no_labelings(); ++i) {
         assert(std::isfinite(msg[i]));
         assert(!std::isnan(msg[i]));
         l[i] += msg[i];
         assert(!std::isnan(l[i]));
         assert(std::isfinite(l[i]));
      }
   }

   template<typename LEFT_FACTOR, typename MSG>
   void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const double omega)
   {
       for(std::size_t i=0; i<l.size(); ++i) {
           assert(std::isfinite(l[i]));
           assert(!std::isnan(l[i]));
       }
      for(std::size_t i=0; i<l.size(); ++i) { assert(!std::isnan(l[i])); }
      msg -= omega*l;
   }

   /*
   template<std::size_t LEFT_IDX, std::size_t... RIGHT_INDICES_REST>
   typename std::enable_if<(LEFT_IDX >= LEFT_LABELINGS::no_labels())>::type
   compute_left_from_right_primal_impl(std::bitset<LEFT_LABELINGS::no_labels()>& l, const std::bitset< RIGHT_LABELINGS::no_labels()>& r)
   {}
   template<std::size_t LEFT_IDX, std::size_t RIGHT_IDX, std::size_t... RIGHT_INDICES_REST>
   typename std::enable_if<(LEFT_IDX < LEFT_LABELINGS::no_labels())>::type
   compute_left_from_right_primal_impl(std::bitset<LEFT_LABELINGS::no_labels()>& l, const std::bitset< RIGHT_LABELINGS::no_labels()>& r)
   {
       l[LEFT_IDX] = r[RIGHT_IDX];
       compute_left_from_right_primal_impl<LEFT_IDX+1, RIGHT_INDICES_REST...>(l,r); 
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
   {
       compute_left_from_right_primal_impl<0, INDICES...>(l.primal(), r.primal()); 
   } 
   */

   template<std::size_t LEFT_IDX, std::size_t... RIGHT_INDICES_REST>
   typename std::enable_if<(LEFT_IDX >= LEFT_LABELINGS::no_labels())>::type
   compute_right_from_left_primal_impl(const std::bitset<LEFT_LABELINGS::no_labels()>& l, std::bitset< RIGHT_LABELINGS::no_labels()>& r)
   {}
   template<std::size_t LEFT_IDX, std::size_t RIGHT_IDX, std::size_t... RIGHT_INDICES_REST>
   typename std::enable_if<(LEFT_IDX < LEFT_LABELINGS::no_labels())>::type
   compute_right_from_left_primal_impl(const std::bitset<LEFT_LABELINGS::no_labels()>& l, std::bitset< RIGHT_LABELINGS::no_labels()>& r)
   {
       r[RIGHT_IDX] = l[LEFT_IDX];
       compute_right_from_left_primal_impl<LEFT_IDX+1, RIGHT_INDICES_REST...>(l,r); 
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
   {
       compute_right_from_left_primal_impl<0, INDICES...>(l.primal(), r.primal());
   }

   template<std::size_t LEFT_IDX, std::size_t... RIGHT_INDICES_REST>
   typename std::enable_if<(LEFT_IDX >= LEFT_LABELINGS::no_labels()),bool>::type
   check_primal_consistency_impl(const std::bitset<LEFT_LABELINGS::no_labels()>& l, const std::bitset< RIGHT_LABELINGS::no_labels()>& r) const
   {
       return true;
   }
   template<std::size_t LEFT_IDX, std::size_t RIGHT_IDX, std::size_t... RIGHT_INDICES_REST>
   typename std::enable_if<(LEFT_IDX < LEFT_LABELINGS::no_labels()),bool>::type
   check_primal_consistency_impl(const std::bitset<LEFT_LABELINGS::no_labels()>& l, const std::bitset< RIGHT_LABELINGS::no_labels()>& r) const
   {
       if( r[RIGHT_IDX] != l[LEFT_IDX] ) { 
           return false; 
       } else {
           return check_primal_consistency_impl<LEFT_IDX+1, RIGHT_INDICES_REST...>(l,r); 
       }
   }


   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
   {
       return check_primal_consistency_impl<0, INDICES...>(l.primal(),r.primal()); 
   }

   template<std::size_t RIGHT_IDX, typename... RIGHT_LABELINGS_REST>
   static typename std::enable_if<(RIGHT_IDX >= RIGHT_LABELINGS::no_labelings())>::type
   print_matching_impl(labelings<RIGHT_LABELINGS_REST...>)
   {}
   template<std::size_t RIGHT_IDX, typename RIGHT_LABELING, typename... RIGHT_LABELINGS_REST>
   static typename std::enable_if<(RIGHT_IDX < RIGHT_LABELINGS::no_labelings())>::type
   print_matching_impl(labelings<RIGHT_LABELING, RIGHT_LABELINGS_REST...>)
   {
      constexpr std::size_t left_label_number = matching_left_labeling<RIGHT_LABELING>();
      std::cout << left_label_number << "," << RIGHT_IDX << "\n";
      print_matching_impl<RIGHT_IDX+1>(labelings<RIGHT_LABELINGS_REST...>{});
   }
   static void print_matching()
   {
      // for each right labeling, print left one that is matched (for debugging purposes)
      print_matching_impl<0>(RIGHT_LABELINGS{});
   }

   template<std::size_t LEFT_LABELING_NO, typename... RIGHT_LABELINGS_REST>
   constexpr static std::size_t no_corresponding_labelings_impl(labelings<RIGHT_LABELINGS_REST...>)
   {
      return 0;
   }
   template<std::size_t LEFT_LABELING_NO, typename RIGHT_LABELING, typename... RIGHT_LABELINGS_REST>
   constexpr static std::size_t no_corresponding_labelings_impl(labelings<RIGHT_LABELING, RIGHT_LABELINGS_REST...>)
   {
      if(matching_left_labeling<RIGHT_LABELING>() == LEFT_LABELING_NO) {
         return 1 + no_corresponding_labelings_impl<LEFT_LABELING_NO>(labelings<RIGHT_LABELINGS_REST...>{});
      } else {
         return no_corresponding_labelings_impl<LEFT_LABELING_NO>(labelings<RIGHT_LABELINGS_REST...>{});
      }
   }
   template<std::size_t LEFT_LABELING_NO>
   constexpr static std::size_t no_corresponding_labelings()
   {
      return no_corresponding_labelings_impl<LEFT_LABELING_NO>(RIGHT_LABELINGS{}); 
   }


   template<std::size_t LEFT_LABELING_NO, std::size_t RIGHT_LABELING_IDX, typename ARRAY_IT, typename... RIGHT_LABELINGS_REST>
   void corresponding_labelings_impl(ARRAY_IT idx, labelings<RIGHT_LABELINGS_REST...>) const
   {
      static_assert(RIGHT_LABELING_IDX == RIGHT_LABELINGS::no_labelings(), "");
      return;
   }
   template<std::size_t LEFT_LABELING_NO, std::size_t RIGHT_LABELING_IDX, typename ARRAY_IT, typename RIGHT_LABELING, typename... RIGHT_LABELINGS_REST>
   void corresponding_labelings_impl(ARRAY_IT it, labelings<RIGHT_LABELING, RIGHT_LABELINGS_REST...>) const
   {
      constexpr std::size_t left_label_number = matching_left_labeling<RIGHT_LABELING>();
      if(left_label_number == LEFT_LABELING_NO) {
         *it = RIGHT_LABELING_IDX;
         corresponding_labelings_impl<LEFT_LABELING_NO, RIGHT_LABELING_IDX+1>(it+1, labelings<RIGHT_LABELINGS_REST...>{});
      } else {
         corresponding_labelings_impl<LEFT_LABELING_NO, RIGHT_LABELING_IDX+1>(it, labelings<RIGHT_LABELINGS_REST...>{});
      } 
   }
   template<std::size_t LEFT_LABELING_NO>
   std::array< std::size_t, no_corresponding_labelings<LEFT_LABELING_NO>() > corresponding_labelings() const
   {
      std::array< std::size_t, no_corresponding_labelings<LEFT_LABELING_NO>() > idx;
      corresponding_labelings_impl<LEFT_LABELING_NO, 0>(idx.begin(), RIGHT_LABELINGS{});
      return idx;
   }


   template<std::size_t LEFT_LABELING_NO, typename EXTERNAL_SOLVER, typename VECTOR>
   typename std::enable_if<(LEFT_LABELING_NO >= LEFT_LABELINGS::no_labelings())>::type
   construct_constraints_impl(EXTERNAL_SOLVER& s, VECTOR left_vars, VECTOR right_vars) const
   {}
   template<std::size_t LEFT_LABELING_NO, typename EXTERNAL_SOLVER, typename VECTOR>
   typename std::enable_if<(LEFT_LABELING_NO < LEFT_LABELINGS::no_labelings())>::type
   construct_constraints_impl(EXTERNAL_SOLVER& s, VECTOR left_vars, VECTOR right_vars) const
   {
      auto right_idx = corresponding_labelings<LEFT_LABELING_NO>();
      std::array<typename EXTERNAL_SOLVER::variable, no_corresponding_labelings<LEFT_LABELING_NO>()> idx;
      for(std::size_t i=0; i<idx.size(); ++i) {
         idx[i] += right_vars[right_idx[i]];
      }
      auto one_active = s.add_at_most_one_constraint(idx.begin(), idx.end());
      s.make_equal(left_vars[LEFT_LABELING_NO], one_active);
      construct_constraints_impl<LEFT_LABELING_NO+1>(s, left_vars, right_vars);
   }
   template<typename EXTERNAL_SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR, typename VECTOR>
   void construct_constraints(EXTERNAL_SOLVER& s, const LEFT_FACTOR& l, VECTOR left_vars, const RIGHT_FACTOR& r, VECTOR right_vars) const
   {
      construct_constraints_impl<0>(s, left_vars, right_vars);
   }

private:
};

} // end namespace LPMP

#endif //  LPMP_LABELING_LIST_FACTOR_HXX

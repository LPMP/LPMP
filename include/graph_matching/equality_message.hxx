#ifndef LPMP_EQUALITY_MESSAGE
#define LPMP_EQUALITY_MESSAGE

#include "config.hxx"
#include <type_traits>

namespace LPMP {

// maximize/minimize to second min/max
// do zrobienia: use breakpoint cost for message updates
class EqualityMessage 
{
public:
   EqualityMessage(const INDEX lV, const INDEX rV)
      : leftVar_(lV), rightVar_(rV)
   {}

   template<typename G1, typename G2>
   void send_message_to_left(const G1& rightPot, G2& msg, const REAL omega)
   {
      MakeFactorUniform(rightPot, msg, rightVar_, omega);
   }
   template<typename G1, typename G2>
   void send_message_to_right(const G1& leftPot, G2& msg, const REAL omega)
   {
      MakeFactorUniform(leftPot, msg, leftVar_, omega);
   }

   template<typename REPAM_ARRAY, typename MSG>
   void MakeFactorUniform(const REPAM_ARRAY& repamPot, MSG& msg, const INDEX var_idx, const REAL omega = 1.0)
   {
      assert(var_idx < repamPot.size());

      // possibly do it differently: search for two second smallest entries and then select first or second one depending upon whether it is rightVar_ or not. Faster?
      REAL min_val = std::numeric_limits<REAL>::max();
      for(INDEX i=0; i<repamPot.size(); ++i) {
         if(i!=var_idx) {
            min_val = std::min(min_val, repamPot[i]);
         }
      }

      msg[0] -= omega*(repamPot[var_idx] - min_val);
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void receive_restricted_message_from_right(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
   {
      if(r.primal() == rightVar_) {
         for(std::size_t i=0; i<l.size(); ++i) {
            if(i != leftVar_) {
               l[i] = std::numeric_limits<REAL>::infinity();
            }
         }
      } else if(r.primal() < r.size()) {
         l[leftVar_] = std::numeric_limits<REAL>::infinity();
      } else {
         //MakeFactorUniform(r, msg, rightVar_);
      }
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void receive_restricted_message_from_left(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
   { 
      if(l.primal() == leftVar_) {
         for(std::size_t i=0; i<r.size(); ++i) {
            if(i != rightVar_) {
               r[i] = std::numeric_limits<REAL>::infinity();
            }
         }
      } else if(l.primal() < l.size()) {
         r[rightVar_] = std::numeric_limits<REAL>::infinity();
      } else {
         //MakeFactorUniform(l, msg, leftVar_);
      } 
   }

   // send all messages of the same type at once
   // possibly templatize this so that case, where all variables of left factor are accessed exatly once, is processed with greater speed (only one minimum searching required then)
   // template specialization is also possible for every dimension accessed except for last one
   // assume that all attached factors have different endpoints

   // do zrobienia: not best encapsulation: we have access to MessageTypeCRTP holding EqualityMessage, not to the EqualityMessage
   // idea: construct proxy object that will give access directly to EqualityMessage, but has no need to construct  the array explicitly
   //
   // code duplication: templatize code

   // for sending multiple messages at once: makes factor uniform by sending all messages at once
   template<typename VAR_ACCESS_OP, typename MSG_ARRAY, typename RIGHT_REPAM>
   static void MakeFactorUniformParallel(VAR_ACCESS_OP var_access_op, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const RIGHT_REPAM& repam, const REAL omega)
   {
      //assert(msgs.size() >= 2); // otherwise calling this method makes no sense, but it can happen for some trivial problems.
      //assert(msgs.size() <= repam.size());
      //assert(msgs.size() == repam.size()-1); // special case of one edge is not assignment in QAP or cosegmentation. For now. Only for hotel and house
      INDEX size = 0;
      for(auto it= msg_begin; it!=msg_end; ++it) ++size;

      assert(omega <= 1.0 + eps);

      // find minimal value of potential over all indices accessed by messages
      REAL min_val_covered = std::numeric_limits<REAL>::max();
      //for(INDEX msg_idx=0; msg_idx<msgs.size(); ++msg_idx) {
      for(auto it= msg_begin; it!=msg_end; ++it) {
         const INDEX var_idx = var_access_op((*it).GetMessageOp());
         //assert(var_idx != repam.size()-1); // this is only valied for assignment problems from house and hotel
         //std::cout << "leftVar = " << leftVar << "\n";
         min_val_covered = std::min(min_val_covered, repam[var_idx]);
      }

      REAL min_val = std::numeric_limits<REAL>::infinity();
      REAL second_min_val = std::numeric_limits<REAL>::infinity();
      for(INDEX i=0; i<repam.size(); ++i) {
         const REAL cur_val = repam[i];
         //std::cout << "cur_val = " << cur_val << "\n";
         if(min_val >= cur_val) { // if two values are equally small, second_min_val should be the lowest value, too
            second_min_val = min_val;
            min_val = cur_val;
         } else if(second_min_val > cur_val) {
            second_min_val = cur_val;
         }
      }
      //assert(std::make_pair(min_val, second_min_val) == SmallestValues<REAL>(repam));

      REAL new_val;  // this value will be taken by the new reparametrized entries
      if(min_val < min_val_covered) { new_val = min_val; }
      else { new_val = second_min_val; }

      //std::cout << "omega_sum = " << omega_sum << ", no messages = " << msgs.size() << ", potential size = " << leftRepam.size() << "\n";
      //std::cout << "min_val_covered = " << min_val_covered << ", min_val = " << min_val << ", second_min_val = " << second_min_val << ", new_val = " << new_val << "\n";
      //const INDEX last_idx = leftRepam.size() - 1;
      //std::cout << "not covered = " << leftRepam[last_idx] << "\n";

      //for(INDEX msg_idx=0; msg_idx<msgs.size(); ++msg_idx, omegaIt++) {
      for(auto it= msg_begin; it!=msg_end; ++it) {
        const INDEX var_idx = var_access_op((*it).GetMessageOp());
        (*it).operator[](0) -= omega*(repam[var_idx] - new_val);
      }
   }

   template<typename RIGHT_REPAM, typename MSG_ARRAY>
   static void SendMessagesToLeft(const RIGHT_REPAM& rightRepam, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
   {
      auto var_access_op = [](const EqualityMessage& msg) -> INDEX { return msg.rightVar_; };
      MakeFactorUniformParallel(var_access_op, msg_begin, msg_end, rightRepam, omega);
   }

   template<typename LEFT_REPAM, typename MSG_ARRAY>
   static void SendMessagesToRight(const LEFT_REPAM& leftRepam, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
   {
      auto var_access_op = [](const EqualityMessage& msg) -> INDEX { return msg.leftVar_; };
      MakeFactorUniformParallel(var_access_op, msg_begin, msg_end, leftRepam, omega);
   }

   template<typename G>
   void RepamLeft(G& leftRepamPot, const REAL msg, const INDEX dim) { 
      assert(dim == 0); 
      leftRepamPot[leftVar_] += msg; 
      assert(!std::isnan(leftRepamPot[leftVar_]));
   }
   template<typename G>
   void RepamRight(G& rightRepamPot, const REAL msg, const INDEX dim) { 
      assert(dim == 0); 
      rightRepamPot[rightVar_] += msg; 
      assert(!std::isnan(rightRepamPot[rightVar_]));
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool
   ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
   {
      if(r.primal() == rightVar_) { 
         const bool ret = (l.primal() != leftVar_);
         l.primal() = leftVar_;
         assert(l.primal() == leftVar_);
         return ret;
      }
      return false;
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool
   ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
   {
      if(l.primal() == leftVar_) { 
         const bool ret = (r.primal() != rightVar_);
         r.primal() = rightVar_;
         assert(r.primal() == rightVar_);
         return ret;
      }
      return false;
   }

   // here it is checked whether labeling on left side and labeling on right side fulfill the constraints of the message
   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
   {
      if(l.primal() == leftVar_) {
         if(r.primal() != rightVar_) std::cout << "kwaskwaskwas1\n";
         return r.primal() == rightVar_;
      }
      if(r.primal() == rightVar_) {
         if(l.primal() != leftVar_) std::cout << "kwaskwaskwas2\n";
         return l.primal() == leftVar_;
      }
      return true;
   }

   template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void construct_constraints(SOLVER& s, LEFT_FACTOR& l, typename SOLVER::vector left_variables, RIGHT_FACTOR& r, typename SOLVER::vector right_variables) const
   {
     s.make_equal(left_variables[leftVar_], right_variables[rightVar_]);
   }



private:
   //do zrobienia: possibly SHORT_INDEX or some 16 bit index (i.e. short unsigned int)
   // certainly 32 bits will be enough.
   const INDEX leftVar_, rightVar_; // variables affected 
};

} // end namespace LPMP

#endif // LPMP_EQUALITY_MESSAGE

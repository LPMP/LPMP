#pragma once

#include <array>
#include "config.hxx"
#include "factors_messages.hxx"
#include "tree_decomposition.hxx"

namespace LPMP {
struct test_factor {
  test_factor(const REAL x, const REAL y)
    : cost(2)
  {
    cost[0] = x;
    cost[1] = y;
  }
  REAL LowerBound() const 
  { 
    return cost.min(); 
  }
  REAL EvaluatePrimal() const 
  { 
    assert(primal < 2);
    return cost[primal]; 
  }

  void MaximizePotentialAndComputePrimal()
  {
    if(primal == std::numeric_limits<INDEX>::max()) {
      primal = std::min_element(cost.begin(),cost.end()) - cost.begin();
    }
    assert(0 <= primal && primal < 2);
  }

  void init_primal() { primal = std::numeric_limits<INDEX>::max(); }

  template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(cost); };
  template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal); }; 

  auto export_variables() { return std::tie(cost); } 

  template<typename ARRAY>
  void apply(ARRAY& a) const
  {
    assert(primal < cost.size());
    a[primal];
  }

  template<typename SOLVER>
  void construct_constraints(SOLVER& s, typename SOLVER::vector v)
  {
    s.add_simplex_constraint(v.begin(), v.end()); 
  }

  template<typename SOLVER>
  void convert_primal(SOLVER& s, typename SOLVER::vector v)
  {
    if(s.solution(v[0])) { primal = 0; }
    else { primal = 1; }
  }

  vector<REAL> cost;
  INDEX primal;
};

struct test_message {

  template<typename LEFT_FACTOR, typename MSG>
  void RepamLeft(LEFT_FACTOR& l, MSG msg)
  {
    assert(msg.size() == 2);
    l.cost[0] += msg[0];
    l.cost[1] += msg[1]; 
  }

  template<typename RIGHT_FACTOR, typename MSG>
  void RepamRight(RIGHT_FACTOR& r, MSG msg)
  {
    assert(msg.size() == 2);
    r.cost[0] += msg[0];
    r.cost[1] += msg[1]; 
  }

  template<typename LEFT_FACTOR, typename MSG>
  void send_message_to_right(const LEFT_FACTOR& l, MSG msg, const REAL omega)
  {
      const auto min = l.cost.min();
      auto m = l.cost;
      for(auto& x : m) { x -= min; }
      msg -= omega*m; 
  }

  template<typename RIGHT_FACTOR, typename MSG>
  void send_message_to_left(const RIGHT_FACTOR& r, MSG msg, const REAL omega)
  {
      const auto min = r.cost.min();
      auto m = r.cost;
      for(auto& x : m) { x -= min; }
      msg -= omega*m; 
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
  {
    r.primal = l.primal;
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
  {
    l.primal = r.primal;
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
  {
    return l.primal == r.primal;
  }


  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, LEFT_FACTOR& l, typename SOLVER::vector v_left, RIGHT_FACTOR& r, typename SOLVER::vector v_right)
  {
    s.make_equal(v_left[0], v_right[0]);
    s.make_equal(v_left[1], v_right[1]);
    //s.make_equal(v_left.begin(), v_left.end(), v_right.begin(), v_right.end()); 
  }
};

struct empty_problem_constructor {
    template<typename SOLVER>
        empty_problem_constructor(SOLVER& s) {}
};

struct test_FMC { // factor message connection
  constexpr static const char* name = "test model";
  using factor = FactorContainer<test_factor, test_FMC, 0>;
  using message = MessageContainer<test_message, 0, 0, message_passing_schedule::left, variableMessageNumber, variableMessageNumber, test_FMC, 0>;
  using FactorList = meta::list<factor>;
  using MessageList = meta::list<message>;
  using problem_constructor = empty_problem_constructor;
};

template<typename LP_TYPE>
void build_test_model(LP_TYPE& lp)
{
  auto* f1 = lp.template add_factor<typename test_FMC::factor>(0,1);
  {
    static_assert( std::is_same_v<typename LP_TYPE::FMC, test_FMC> );
    factor_tree<test_FMC> t1;
    auto* f2 = lp.template add_factor<typename test_FMC::factor>(1,0);
    auto* f3 = lp.template add_factor<typename test_FMC::factor>(0,0);
    auto* m12 = lp.template add_message<typename test_FMC::message>(f1,f2);
    auto* m13 = lp.template add_message<typename test_FMC::message>(f1,f3);
    t1.add_message(m12, Chirality::left);
    t1.add_message(m13, Chirality::left);
    lp.add_tree(t1);
  }

  {
    factor_tree<test_FMC> t2;
    auto* f2 = lp.template add_factor<typename test_FMC::factor>(1,0);
    auto* f3 = lp.template add_factor<typename test_FMC::factor>(0,0);
    auto* m12 = lp.template add_message<typename test_FMC::message>(f1,f2);
    auto* m23 = lp.template add_message<typename test_FMC::message>(f2,f3);
    t2.add_message(m12, Chirality::right);
    t2.add_message(m23, Chirality::left);
    lp.add_tree(t2);
  }

  {
    factor_tree<test_FMC> t3;
    auto* f2 = lp.template add_factor<typename test_FMC::factor>(1,0);
    auto* f3 = lp.template add_factor<typename test_FMC::factor>(0,0);
    auto* m12 = lp.template add_message<typename test_FMC::message>(f1,f2);
    auto* m23 = lp.template add_message<typename test_FMC::message>(f2,f3);
    t3.add_message(m12, Chirality::right);
    t3.add_message(m23, Chirality::left);
    lp.add_tree(t3);
  }
}

} // namespace LPMP 

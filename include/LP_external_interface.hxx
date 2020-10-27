#ifndef LPMP_LP_EXTERNAL_INTERFACE_HXX
#define LPMP_LP_EXTERNAL_INTERFACE_HXX

#include "DD_ILP.hxx"
#include "LP.h"
#include "external_solver_interface.hxx"

// interface to DD_ILP object which builds up the given LPMP problem for various other solvers.
// there are two functions a factor must provide, so that the export can take place:
//   template<typename EXTERNAL_SOLVER>
//   void construct_constraints(EXTERNAL_SOLVER& s, TYPES...) const
// and
//   template<typename EXTERNAL_SOLVER>
//   void convert_to_primal(EXTERNAL_SOLVER& s, TYPES...)
// where TYPES is a list of parameters of type EXTERNAL_SOLVER::{variable|vector|matrix|tensor} whose number and type matches the respective arguments given to export_variables()
// likewise, a message must also provide a function
// template<typename EXTERNAL_SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
// void construct_constraints(EXTERNAL_SOLVER&s, LEFT_FACTOR& left, LEFT_TYPES..., RIGHT_FACTOR& right, RIGHT_TYPES...)

namespace LPMP {

template<typename BASE_EXTERNAL_SOLVER, typename BASE_LP_SOLVER>
class LP_external_solver : public BASE_LP_SOLVER {
public:
  using external_solver = DD_ILP::external_solver_interface<BASE_EXTERNAL_SOLVER>;

  using BASE_LP_SOLVER::BASE_LP_SOLVER;

  void solve()
  {
    construct();
    if (!s_.solve())
      throw std::runtime_error("External solver was unable to optimize model.");

    s_.init_variable_loading();
    for (auto* f : this->f_) {
      f->convert_primal(s_);
    }

    assert(this->CheckPrimalConsistency());
  }

  void write_to_file(const std::string& filename) {
    construct();
    s_.write_to_file(filename);
  }

  const external_solver& get_external_solver() const { return s_; }

private:
  void construct() {
    s_.init_variable_loading();

    if (external_variable_counter_.size() > 0)
      return; // already constructed

    // Can't use `for_each_factor` here, as the order is different than `f_`
    // and we rely on `factor_address_to_index_` which gets updated in
    // `add_factor`.
    for (auto* f : this->f_) {
      external_variable_counter_.push_back(s_.get_variable_counters());
      f->construct_constraints(s_);
    }

    this->for_each_message([&](auto* m) {
      const INDEX left_factor_no = this->factor_address_to_index_[m->GetLeftFactor()];
      assert(left_factor_no < this->number_of_factors() && left_factor_no < external_variable_counter_.size());

      const INDEX right_factor_no = this->factor_address_to_index_[m->GetRightFactor()];
      assert(right_factor_no < this->number_of_factors() && right_factor_no < external_variable_counter_.size());

      m->construct_constraints(s_, external_variable_counter_[left_factor_no], external_variable_counter_[right_factor_no]);
    });

    s_.init_variable_loading();
    for (auto* f : this->f_) // Some note regarding `f_` as above.
      f->load_costs(s_);
  }

  external_solver s_;

  std::vector<typename DD_ILP::variable_counters> external_variable_counter_;
};


} // end namespace LPMP

#endif // LPMP_LP_EXTERNAL_INTERFACE_HXX

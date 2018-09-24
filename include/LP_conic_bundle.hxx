#ifndef LPMP_LP_CONIC_BUNDLE_HXX
#define LPMP_LP_CONIC_BUNDLE_HXX

#include "tree_decomposition.hxx"
#include "CBSolver.hxx"

namespace LPMP {

template<typename FMC_TYPE>
class LP_conic_bundle : public LP_with_trees<FMC_TYPE, Lagrangean_factor_star, LP_conic_bundle<FMC_TYPE> >, public ConicBundle::FunctionOracle {
public:
   using FMC = FMC_TYPE;
   using LP_with_trees<FMC, Lagrangean_factor_star, LP_conic_bundle<FMC> >::LP_with_trees;

   void construct_decomposition()
   {
      std::cout << "no of Lagrangean vars: " << this->no_Lagrangean_vars() << "\n";

      //cb_solver_.set_max_bundlesize(*this,10);
      //cb_solver_.set_eval_limit(1250); // recommendation by Jan Kuske

      // set up problem for conic bundle
      cb_solver_.init_problem(this->no_Lagrangean_vars());
      cb_solver_.set_out(&std::cout,1);
      cb_solver_.add_function(*this);
      cb_solver_.do_descent_step(); // we need to setup the solver like this as well
   }

   void optimize_decomposition()
   {
      cb_solver_.do_descent_step();

      if(cb_solver_.termination_code()) {
         std::cout << "\nconic bundle has terminated\n\n";
         cb_solver_.print_termination_code(std::cout);
      }
   }

   int evaluate(const ConicBundle::DVector& x, double relprec, 
         double& objective_value,
         ConicBundle::DVector&  cut_vals,std::vector<ConicBundle::DVector>&  subgradients,
         std::vector<ConicBundle::PrimalData*>&  primal_solutions,
         ConicBundle::PrimalExtender*&)
   {
      // load Lagrangean variables
      this->add_weights(&x[0], -1.0);

      // compute subgradient
      objective_value = 0.0;
      ConicBundle::DVector subg(x.size(), 0.0); // this is not so nice!
      for(auto& t : this->trees_) {
         t.solve();
         t.compute_mapped_subgradient(subg);
         objective_value -= t.primal_cost();
      }
      cut_vals.push_back(objective_value);
      subgradients.push_back(subg);

      // remove Lagrangean variables again
      this->add_weights(&x[0], +1.0);

      return 0;
   }

   REAL decomposition_lower_bound() 
   {
      assert(cb_solver_.get_dim() == this->no_Lagrangean_vars() && cb_solver_.get_dim() > 0);
      ConicBundle::DVector Lagrangean_vars;
      cb_solver_.get_center(Lagrangean_vars);
      // load Lagrangean variables
      this->add_weights(&Lagrangean_vars[0], -1.0);

      const REAL lb = LP_with_trees<FMC, Lagrangean_factor_star, LP_conic_bundle<FMC> >::decomposition_lower_bound();
      // remove Lagrangean variables again
      this->add_weights(&Lagrangean_vars[0], +1.0);

      return lb;
   }

private:
   ConicBundle::CBSolver cb_solver_;

};

} // end namespace LPMP

#endif // LPMP_LP_CONIC_BUNDLE_HXX

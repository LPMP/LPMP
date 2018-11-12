#ifndef LPMP_FACTOR_CONTAINER_INTERFACE
#define LPMP_FACTOR_CONTAINER_INTERFACE

#include "message_passing_schedule.hxx"
#include "two_dimensional_variable_array.hxx"
#include "serialization.hxx"
#include "DD_ILP.hxx"
#include "config.hxx"

namespace LPMP {

using weight_array = two_dim_variable_array<REAL>;
using weight_slice = two_dim_variable_array<REAL>::ArrayAccessObject;
using receive_array = two_dim_variable_array<unsigned char>;
using receive_slice = two_dim_variable_array<unsigned char>::ArrayAccessObject;

// pure virtual base class for factor container used by LP class
class FactorTypeAdapter
{
public:
   virtual ~FactorTypeAdapter() {}
   virtual FactorTypeAdapter* clone() const = 0;
   virtual void update_factor_uniform(const REAL leave_weight) = 0;
   virtual void UpdateFactor(const weight_slice omega, const receive_slice receive_mask) = 0;
   virtual void update_factor_residual(const weight_slice omega, const receive_slice receive_mask) = 0;
   virtual void UpdateFactorPrimal(const weight_slice& omega, const receive_slice& receive_mask, const INDEX iteration) = 0;
   // TODO: remove below two functions. Can be handled via message passing schedule
   virtual bool SendsMessage(const INDEX msg_idx) const = 0;
   virtual bool ReceivesMessage(const INDEX msg_idx) const = 0;
   virtual bool FactorUpdated() const = 0; // does calling UpdateFactor do anything? If no, it need not be called while in ComputePass, saving time.
   virtual std::size_t no_messages() const = 0;
   virtual std::size_t no_send_messages() const = 0;
   virtual std::size_t no_receive_messages() const = 0;
   virtual REAL LowerBound() const = 0;
   virtual void init_primal() = 0;
   virtual void MaximizePotentialAndComputePrimal() = 0;
   virtual void propagate_primal_through_messages() = 0;
   virtual bool check_primal_consistency() = 0;

   // for use in tree decomposition:
   // for writing primal solution into subgradient
   // return value is size of subgradient
   virtual INDEX subgradient(double* w, const REAL sign) = 0;
   virtual REAL dot_product(double* w) = 0;

   // for reading reparametrization/labeling out of factor
   virtual void serialize_dual(save_archive&) = 0;
   virtual void serialize_primal(save_archive&) = 0;
   // for writing reparametrization/labeling into factor
   virtual void serialize_dual(load_archive&) = 0;
   virtual void serialize_primal(load_archive&) = 0;
   // for determining size of archive
   virtual void serialize_dual(allocate_archive&) = 0;
   virtual void serialize_primal(allocate_archive&) = 0;
   // for adding weights in Frank Wolfe algorithm
   virtual void serialize_dual(addition_archive&) = 0;

   virtual void divide(const REAL val) = 0; // divide potential by value
   virtual void set_to_value(const REAL val) = 0; // set potential to given value
   virtual void add(FactorTypeAdapter*) = 0; // add potential values of other factor

   virtual INDEX dual_size() = 0;
   virtual INDEX dual_size_in_bytes() = 0;
   virtual INDEX primal_size_in_bytes() = 0;

   // do zrobienia: this function is not needed. Evaluation can be performed automatically
   virtual REAL EvaluatePrimal() const = 0;

   // external ILP-interface
   virtual void construct_constraints(DD_ILP::external_solver_interface<DD_ILP::problem_export>& solver) = 0;
   virtual void load_costs(DD_ILP::external_solver_interface<DD_ILP::problem_export>& solver) = 0; 
   virtual void convert_primal(DD_ILP::external_solver_interface<DD_ILP::problem_export>& solver) = 0;

   virtual void construct_constraints(DD_ILP::external_solver_interface<DD_ILP::sat_solver>& solver) = 0;
   virtual void load_costs(DD_ILP::external_solver_interface<DD_ILP::sat_solver>& solver) = 0;
   virtual void convert_primal(DD_ILP::external_solver_interface<DD_ILP::sat_solver>& solver) = 0; 

#ifdef DD_ILP_WITH_GUROBI
   virtual void construct_constraints(DD_ILP::external_solver_interface<DD_ILP::gurobi_interface>& solver) = 0;
   virtual void load_costs(DD_ILP::external_solver_interface<DD_ILP::gurobi_interface>& solver) = 0;
   virtual void convert_primal(DD_ILP::external_solver_interface<DD_ILP::gurobi_interface>& solver) = 0;
#endif

   // estimate of how long a factor update will take
   virtual INDEX runtime_estimate() = 0;

   virtual std::vector<FactorTypeAdapter*> get_adjacent_factors() const = 0;

   struct message_trait {
       FactorTypeAdapter* adjacent_factor;
       Chirality chirality; // is f on the left or on the right
       message_passing_schedule_factor_view::type mps;
   };
   virtual std::vector<message_trait> get_messages() const = 0;
};

} // namespace LPMP

#endif // LPMP_FACTOR_CONTAINER_INTERFACE

#pragma once

#include <vector>
#include <cassert>
#include <future>
#include "cut_base/cut_base_triplet_constructor.hxx"
#include "max_cut_instance.hxx"
#include "max_cut_cycle_packing.h"
#include "max_cut_local_search.h"
#include "graph.hxx"
#include "max_cut_greedy_additive_edge_contraction.h"
#include "max_cut_sahni_gonzalez.h"
#include "mrf/binary_MRF_instance.hxx"
#include "mrf/qpbo_factor.hxx"
#include "max_cut/transform_binary_MRF.hxx"

namespace LPMP {

template<
   class FACTOR_MESSAGE_CONNECTION, 
   typename UNARY_FACTOR, typename TRIPLET_FACTOR,
   typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2
   >
class max_cut_triplet_constructor : public cut_base_triplet_constructor<
                                    FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2> 
{
public:
   using type = max_cut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>;
   using FMC = FACTOR_MESSAGE_CONNECTION;
   using base_constructor = cut_base_triplet_constructor<
       FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>;

   using weighted_graph = graph<double>;

   template<typename SOLVER>
       max_cut_triplet_constructor(SOLVER& s);

   std::size_t find_violated_cycles(const std::size_t max_triplets_to_add);

   void construct(const max_cut_instance& mc);
   void construct(const binary_MRF_instance& instance);

   bool CheckPrimalConsistency() const;
   std::size_t Tighten(const std::size_t no_constraints);
   //static std::vector<char> round(std::vector<typename base_constructor::edge> edges);
   static max_cut_edge_labeling round(const max_cut_instance input, const std::string method, const max_cut_instance& original_model);
   void ComputePrimal();
   void Begin();
   void End();
private:
    std::future<max_cut_edge_labeling> primal_result_handle_;
    std::future<cycle_packing> cycle_inequalities_tighten_handle_;
    max_cut_instance original_model;

    TCLAP::ValueArg<std::string> rounding_method_arg_;
};

template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
    template<typename SOLVER>
max_cut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::max_cut_triplet_constructor(SOLVER& s)
    :
        rounding_method_arg_("", "maxCutRounding", "method for rounding primal solution", false, "gaec", "{gaec|sahni_gonzalez}", s.get_cmd()),
        base_constructor(s)
{}

template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
   void max_cut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::construct(const max_cut_instance& mc)
   {
      for(const auto& e : mc.edges())
         this->add_edge_factor(e[0], e[1], e.cost);
      this->no_original_edges_ = this->unary_factors_vector_.size();
      this->lp_->add_to_constant(mc.constant());
   }

    template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
void max_cut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::construct(const binary_MRF_instance& instance)
{
    std::cout << "original lower bound = " << instance.lower_bound() << "\n";
    qpbo_factor q(instance);
    const binary_MRF_instance reduced_instance = q.get_reduced_problem();
    std::cout << "lb of reduced mrf = " << reduced_instance.lower_bound() << "\n";
    const binary_Potts_instance reduced_Potts_instance = transform_binary_MRF_to_Potts(reduced_instance);
    const max_cut_instance reduced_max_cut = transform_binary_Potts_to_max_cut(reduced_Potts_instance);
    if(diagnostics()) {
        std::cout << "Solved QPBO, lower bound = " << q.LowerBound() << ", #original variables = " << instance.unaries.size() << ", #reduced variables = " << reduced_instance.unaries.size() << "\n";
    }
    this->construct(reduced_max_cut); 
    this->lp_->add_to_constant(q.LowerBound());
}

template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
   bool max_cut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::CheckPrimalConsistency() const
   {
       const max_cut_edge_labeling labeling = this->template export_edge_labeling<max_cut_edge_labeling>();
       const auto model = this->template export_edges<max_cut_instance>();
       return labeling.check_primal_consistency(model); 
   }

    template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
    std::size_t max_cut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::Tighten(const std::size_t no_constraints)
{
    /*
    std::size_t no_cycles_added = 0;
    if(cycle_inequalities_tighten_handle_.valid() && cycle_inequalities_tighten_handle_.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
        if(debug())
            std::cout << "add found cycle inequalities\n";
        cycle_packing cp = cycle_inequalities_tighten_handle_.get();
        for(std::size_t c=0; c<cp.no_cycles(); ++c) {
            const auto [cycle_begin, cycle_end] = cp.get_cycle(c);
            const double cycle_weight = cp.get_cycle_weight(c);
            this->triangulate_cycle(cycle_begin, cycle_end, cycle_weight);
        } 
        if(debug())
            std::cout << "Added " << cp.no_cycles() << " cycles through cycle packing\n";
        no_cycles_added = cp.no_cycles(); 
    }

    if(!cycle_inequalities_tighten_handle_.valid() || cycle_inequalities_tighten_handle_.wait_for(std::chrono::seconds(0)) == std::future_status::deferred) {
        if(debug())
            std::cout << "search for violated cycle inequalities\n";
        this->send_messages_to_edges(); // TODO: check if this is always helpful or only for triplet tightening
        const auto mc = this->template export_edges<max_cut_instance>();
        cycle_inequalities_tighten_handle_ = std::async(std::launch::async, compute_max_cut_cycle_packing, std::move(mc));
    }

    return no_cycles_added;
    */

    this->send_messages_to_edges();
    const auto mc = this->template export_edges<max_cut_instance>();
    cycle_packing cp = compute_max_cut_cycle_packing(mc);
    for(std::size_t c=0; c<cp.no_cycles(); ++c) {
        const auto [cycle_begin, cycle_end] = cp.get_cycle(c);
        const double cycle_weight = cp.get_cycle_weight(c);
        this->triangulate_cycle(cycle_begin, cycle_end, cycle_weight);
    } 
    if(debug())
        std::cout << "Added " << cp.no_cycles() << " cycles through cycle packing\n";
    return cp.no_cycles();
}

    template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
    max_cut_edge_labeling max_cut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::round(const max_cut_instance input, const std::string method, const max_cut_instance& original_model)
{
    if(method == "gaec") {
        const max_cut_edge_labeling sol = greedy_additive_edge_contraction(input);
        if(original_model.no_nodes() > 0) {
            std::cout << "local search postprocessing\n";
            const auto node_sol = sol.transform_to_node_labeling(input);
            max_cut_local_search ls(original_model, node_sol);
            ls.perform_swaps();
            const max_cut_node_labeling improved_sol = ls.get_labeling();
            return max_cut_edge_labeling(input,improved_sol);
        } else {
            return sol;
        }
    } else if(method == "sahni_gonzalez") {
        max_cut_node_labeling sol = max_cut_sahni_gonzalez_3(input);
        max_cut_local_search ls(input, sol);
        ls.perform_swaps();
        sol = ls.get_labeling();
        if(original_model.no_nodes() > 0) {
            std::cout << "local search postprocessing\n";
            max_cut_local_search ls(original_model, sol);
            ls.perform_swaps();
            const max_cut_node_labeling improved_sol = ls.get_labeling();
            return max_cut_edge_labeling(input,improved_sol);
        } else {
            return max_cut_edge_labeling(input,sol);
        }
    } else {
        throw std::runtime_error(std::string("rounding method not recognized: ") + method);
    }
}

    template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
    void max_cut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::ComputePrimal()
{
    this->send_messages_to_edges(); // If we do it only when needed, we will get varying results depending on how fast primal rounding is.

    if(primal_result_handle_.valid() && primal_result_handle_.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
        if(debug()) 
            std::cout << "read in primal max_cut solution\n";
        auto max_cut_sol = primal_result_handle_.get();
        this->write_labeling_into_factors(max_cut_sol); 
    }

    if(!primal_result_handle_.valid() || primal_result_handle_.wait_for(std::chrono::seconds(0)) == std::future_status::deferred) {
        if(debug()) 
            std::cout << "export max_cut problem for rounding\n";

        const auto instance = this->template export_edges<max_cut_instance>();
        primal_result_handle_ = std::async(std::launch::async, round, std::move(instance), rounding_method_arg_.getValue(), original_model);
    } 
}

    template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
void max_cut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::Begin()
{
    if(this->number_of_triplets() == 0) { // get original model and perform local search on it to improve rounded solution
        if(debug())
            std::cout << "storing original model for local search postprocessing\n";
        for(const auto& e : this->unary_factors_vector_)
            original_model.add_edge(e.first[0], e.first[1], (*e.second->get_factor())[0]);
    } else {
        if(debug())
            std::cout << "model has triplets, no local search postprocessing possible\n";
    }

    base_constructor::Begin();
    const auto instance = this->template export_edges<max_cut_instance>();
    primal_result_handle_ = std::async(std::launch::async, round, instance, rounding_method_arg_.getValue(), original_model);

    // compute cycle packing and add returned cycles to problem formulation. Also reparametrize edges (possibly do not do this?)
    cycle_packing cp = compute_max_cut_cycle_packing(instance);
    for(std::size_t c=0; c<cp.no_cycles(); ++c) {
        const auto [cycle_begin, cycle_end] = cp.get_cycle(c);
        const double cycle_weight = cp.get_cycle_weight(c);
        this->triangulate_cycle(cycle_begin, cycle_end, cycle_weight);
    } 
}

    template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
void max_cut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::End()
{
    base_constructor::End();
    if(primal_result_handle_.valid()) {
        primal_result_handle_.wait();
        if(debug()) 
            std::cout << "read in primal max_cut solution\n";
        auto max_cut_sol = primal_result_handle_.get();
        this->write_labeling_into_factors(max_cut_sol); 
    }
}

} // namespace LPMP 

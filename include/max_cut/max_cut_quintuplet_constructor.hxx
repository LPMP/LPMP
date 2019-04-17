#pragma once

#include "cut_base/cut_base_quintuplet_constructor.hxx"
#include "max_cut_odd_bicycle_wheel_packing.h"
#include <future>

namespace LPMP {

    template<typename TRIPLET_CONSTRUCTOR, 
        typename QUINTUPLET_FACTOR,
        typename... MSGS
            >
            class max_cut_quintuplet_constructor : public cut_base_triplet_quintuplet_constructor<TRIPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, MSGS...>
            {
                public:
                    using base_constructor = cut_base_triplet_quintuplet_constructor<TRIPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, MSGS...>;

                    using quintuplet_factor_container = QUINTUPLET_FACTOR;

                    template<typename SOLVER>
                        max_cut_quintuplet_constructor(SOLVER& s) : base_constructor(s) {}

                    std::size_t Tighten(const std::size_t no_constraints_to_add);
                    void ComputePrimal();

                private:
                    std::future<odd_bicycle_wheel_packing> odd_bicycle_wheel_inequalities_tighten_handle_;
            };

    template<typename TRIPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename... MSGS>
        std::size_t max_cut_quintuplet_constructor<TRIPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, MSGS...>::Tighten(const std::size_t no_constraints_to_add)
        {
            /*
            std::size_t no_quintuplets_added = 0;
            if(odd_bicycle_wheel_inequalities_tighten_handle_.valid() && odd_bicycle_wheel_inequalities_tighten_handle_.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
                const odd_bicycle_wheel_packing obwp = odd_bicycle_wheel_inequalities_tighten_handle_.get();

                if(debug())
                    std::cout << "found " << obwp.no_odd_bicycle_wheels() << " odd bicycle wheels\n";
                const std::size_t no_quintuplets_before = this->quintuplet_factors().size();
                for(std::size_t c=0; c<obwp.no_odd_bicycle_wheels(); ++c) {
                    const auto [cycle_begin, cycle_end] = obwp.get_cycle(c);
                    const double odd_bicycle_wheel_weight = obwp.get_odd_bicycle_wheel_weight(c);
                    const auto axle = obwp.get_axle(c);
                    this->triangulate_odd_bicycle_wheel(axle, cycle_begin, cycle_end, odd_bicycle_wheel_weight);
                }
                const std::size_t no_quintuplets_after = this->quintuplet_factors().size();
                if(debug())
                    std::cout << "Added " << no_quintuplets_after - no_quintuplets_before << " quintuplets\n";

                no_quintuplets_added = base_constructor::Tighten(no_constraints_to_add);
            }

            if(!odd_bicycle_wheel_inequalities_tighten_handle_.valid() || odd_bicycle_wheel_inequalities_tighten_handle_.wait_for(std::chrono::seconds(0)) == std::future_status::deferred) {
                this->send_messages_to_triplets();
                triplet_max_cut_instance mc = this->template export_triplets<triplet_max_cut_instance>();
                odd_bicycle_wheel_inequalities_tighten_handle_ = std::async(std::launch::async, compute_max_cut_odd_bicycle_wheel_packing, std::move(mc));
            }

            return no_quintuplets_added;
            */

            this->send_messages_to_triplets();
            triplet_max_cut_instance mc = this->template export_triplets<triplet_max_cut_instance>();
            const odd_bicycle_wheel_packing obwp = compute_max_cut_odd_bicycle_wheel_packing(mc);
            if(debug())
                std::cout << "found " << obwp.no_odd_bicycle_wheels() << " odd bicycle wheels\n";
            const std::size_t no_quintuplets_before = this->quintuplet_factors().size();
            for(std::size_t c=0; c<obwp.no_odd_bicycle_wheels(); ++c) {
                const auto [cycle_begin, cycle_end] = obwp.get_cycle(c);
                const double odd_bicycle_wheel_weight = obwp.get_odd_bicycle_wheel_weight(c);
                const auto axle = obwp.get_axle(c);
                this->triangulate_odd_bicycle_wheel(axle, cycle_begin, cycle_end, odd_bicycle_wheel_weight);
            }
            const std::size_t no_quintuplets_after = this->quintuplet_factors().size();
            if(debug())
                std::cout << "Added " << no_quintuplets_after - no_quintuplets_before << " quintuplets\n";

            const std::size_t base_constraints_added = base_constructor::Tighten(no_constraints_to_add);
            
            return base_constraints_added + obwp.no_odd_bicycle_wheels();
        }

    template<typename TRIPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename... MSGS>
        void max_cut_quintuplet_constructor<TRIPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, MSGS...>::ComputePrimal()
        {
            //this->send_messages_to_triplets();
            base_constructor::ComputePrimal(); 
        }

} // namespace LPMP

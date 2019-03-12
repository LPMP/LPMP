#pragma once

#include "cut_base/cut_base_quadruplet_constructor.hxx"
#include "multicut_odd_wheel_packing.h"

namespace LPMP {

    template<typename BASE_CONSTRUCTOR,
        typename QUADRUPLET_FACTOR,
        typename TRIPLET_QUADRUPLET_MESSAGE_012,
        typename TRIPLET_QUADRUPLET_MESSAGE_013,
        typename TRIPLET_QUADRUPLET_MESSAGE_023,
        typename TRIPLET_QUADRUPLET_MESSAGE_123>
            class multicut_quadruplet_constructor : public cut_base_quadruplet_constructor<BASE_CONSTRUCTOR, QUADRUPLET_FACTOR, TRIPLET_QUADRUPLET_MESSAGE_012, TRIPLET_QUADRUPLET_MESSAGE_013, TRIPLET_QUADRUPLET_MESSAGE_023, TRIPLET_QUADRUPLET_MESSAGE_123>
        {
            public:
                using FMC = typename BASE_CONSTRUCTOR::FMC;

                using base_constructor = cut_base_quadruplet_constructor<BASE_CONSTRUCTOR, QUADRUPLET_FACTOR, TRIPLET_QUADRUPLET_MESSAGE_012, TRIPLET_QUADRUPLET_MESSAGE_013, TRIPLET_QUADRUPLET_MESSAGE_023, TRIPLET_QUADRUPLET_MESSAGE_123>;

                using base_constructor::base_constructor;

                /*
                   quadruplet labeling order:
                   1,1,0,1,0,0
                   1,0,1,0,1,0
                   0,1,1,0,0,1
                   0,0,0,1,1,1
                   0,1,1,1,1,0
                   1,0,1,1,0,1
                   1,1,0,0,1,1
                   0,1,1,1,1,1
                   1,0,1,1,1,1
                   1,1,0,1,1,1
                   1,1,1,0,1,1
                   1,1,1,1,0,1
                   1,1,1,1,1,0
                   1,1,1,1,1,1
                   */ 
                template<typename ITERATOR>
                    QUADRUPLET_FACTOR* add_higher_order_quadruplet( const std::size_t i0, const std::size_t i1, const std::size_t i2, const std::size_t i3, ITERATOR cost_begin, ITERATOR cost_end);

                std::size_t Tighten(const std::size_t no_constraints_to_add);
                void ComputePrimal();
        };


    // add to cut base constructor
    template<typename BASE_CONSTRUCTOR, typename QUADRUPLET_FACTOR, typename TRIPLET_QUADRUPLET_MESSAGE_012, typename TRIPLET_QUADRUPLET_MESSAGE_013, typename TRIPLET_QUADRUPLET_MESSAGE_023, typename TRIPLET_QUADRUPLET_MESSAGE_123>
        template<typename ITERATOR>
        QUADRUPLET_FACTOR* multicut_quadruplet_constructor<BASE_CONSTRUCTOR, QUADRUPLET_FACTOR, TRIPLET_QUADRUPLET_MESSAGE_012, TRIPLET_QUADRUPLET_MESSAGE_013, TRIPLET_QUADRUPLET_MESSAGE_023, TRIPLET_QUADRUPLET_MESSAGE_123>::add_higher_order_quadruplet
        ( const std::size_t i0, const std::size_t i1, const std::size_t i2, const std::size_t i3,
          ITERATOR cost_begin, ITERATOR cost_end)
        {
            auto* f = this->add_odd_3_wheel_factor(i0,i1,i2,i3);
            assert(std::distance(cost_begin, cost_end) == 15);
            this->lp_->add_to_constant(*cost_begin);
            auto* t = f->get_factor();

            std::size_t i = 0;
            for(auto it = cost_begin+1; it != cost_end; ++it, ++i) {
                (*t)[i] = *it - *cost_begin; 
            }

            return f;
        }

    template<typename BASE_CONSTRUCTOR, typename QUADRUPLET_FACTOR, typename TRIPLET_QUADRUPLET_MESSAGE_012, typename TRIPLET_QUADRUPLET_MESSAGE_013, typename TRIPLET_QUADRUPLET_MESSAGE_023, typename TRIPLET_QUADRUPLET_MESSAGE_123>
        std::size_t multicut_quadruplet_constructor<BASE_CONSTRUCTOR, QUADRUPLET_FACTOR, TRIPLET_QUADRUPLET_MESSAGE_012, TRIPLET_QUADRUPLET_MESSAGE_013, TRIPLET_QUADRUPLET_MESSAGE_023, TRIPLET_QUADRUPLET_MESSAGE_123>::Tighten(const std::size_t no_constraints_to_add)
        {
            // TODO: search asynchronously
            this->send_messages_to_triplets();
            triplet_multicut_instance mc = this->template export_triplets<triplet_multicut_instance>();
            auto owp = compute_multicut_odd_wheel_packing(mc);
            if(debug())
                std::cout << "found " << owp.no_odd_wheels() << " odd wheels\n";
            for(std::size_t c=0; c<owp.no_odd_wheels(); ++c) {
                const auto [cycle_begin, cycle_end] = owp.get_cycle(c);
                const double odd_wheel_weight = owp.get_odd_wheel_weight(c);
                const std::size_t center_node = owp.get_center_node(c);
                this->triangulate_odd_wheel(center_node, cycle_begin, cycle_end, odd_wheel_weight);
            } 

            const std::size_t base_constraints_added = BASE_CONSTRUCTOR::Tighten(no_constraints_to_add);

            return base_constraints_added + owp.no_odd_wheels();
        }

    template<typename BASE_CONSTRUCTOR, typename QUADRUPLET_FACTOR, typename TRIPLET_QUADRUPLET_MESSAGE_012, typename TRIPLET_QUADRUPLET_MESSAGE_013, typename TRIPLET_QUADRUPLET_MESSAGE_023, typename TRIPLET_QUADRUPLET_MESSAGE_123>
        void multicut_quadruplet_constructor<BASE_CONSTRUCTOR, QUADRUPLET_FACTOR, TRIPLET_QUADRUPLET_MESSAGE_012, TRIPLET_QUADRUPLET_MESSAGE_013, TRIPLET_QUADRUPLET_MESSAGE_023, TRIPLET_QUADRUPLET_MESSAGE_123>::ComputePrimal()
        {
            this->send_messages_to_triplets();
            base_constructor::ComputePrimal();
        }

} // namespace LPMP

#ifndef LPMP_FACTOR_PARTITIONS_STORAGE_HXX
#define LPMP_FACTOR_PARTITIONS_STORAGE_HXX

#include <vector>
#include <array>
#include "two_dimensional_variable_array.hxx"
#include "factor_container_interface.h"
#include "union_find.hxx"

namespace LPMP {

class factors_partition_storage {
public:
   // methods for staged optimization
   void put_in_same_partition(FactorTypeAdapter* f1, FactorTypeAdapter* f2) { factor_partition_valid_ = false; partition_graph.push_back({f1,f2}); }
   void

   void construct_factor_partition();
   void construct_overlapping_factor_partition();

   template<typename PARTITION_ITERATOR, typename INTRA_PARTITION_FACTOR_ITERATOR>
   void construct_forward_pushing_weights(PARTITION_ITERATOR partition_begin, PARTITION_ITERATOR partition_end, std::vector<weight_array>& omega_partition, std::vector<receive_array>& receive_mask_partition, INTRA_PARTITION_FACTOR_ITERATOR factor_iterator_getter);

   void compute_partition_pass(const std::size_t no_passes);
   void compute_overlapping_partition_pass(const std::size_t no_passes);

private:
   template<typename ITERATOR_1, typename ITERATOR_2>
   std::vector<FactorTypeAdapter*> concatenate_factors(ITERATOR_1 f1_begin, ITERATOR_1 f1_end, ITERATOR_2 f2_begin, ITERATOR_2 f2_end);

   // for staged optimization: partition factors that are updated and run multiple rounds of optimization on each component of the partition followed by pushing messages to the next component.
   std::vector<std::array<FactorTypeAdapter*,2>> partition_graph;
   bool factor_partition_valid_ = false;
   two_dim_variable_array<FactorTypeAdapter*> factor_partition_; // already sorted
   std::vector<weight_array> omega_partition_forward_, omega_partition_backward_;
   std::vector<receive_array> receive_mask_partition_forward_, receive_mask_partition_backward_;

   std::vector<weight_array> omega_partition_forward_pass_push_forward_, omega_partition_backward_pass_push_forward_;
   std::vector<weight_array> omega_partition_forward_pass_push_backward_, omega_partition_backward_pass_push_backward_;
   std::vector<receive_array> receive_mask_partition_forward_pass_push_forward_, receive_mask_partition_backward_pass_push_forward_;
   std::vector<receive_array> receive_mask_partition_forward_pass_push_backward_, receive_mask_partition_backward_pass_push_backward_;

   std::vector<weight_array> omega_partition_forward_pass_push_;
   std::vector<weight_array> omega_partition_backward_pass_push_;
   std::vector<receive_array> receive_mask_partition_forward_pass_push_;
   std::vector<receive_array> receive_mask_partition_backward_pass_push_;

   bool overlapping_factor_partition_valid_ = false;
   std::vector<weight_array> omega_overlapping_partition_forward_;
   std::vector<weight_array> omega_overlapping_partition_backward_;
   std::vector<receive_array> receive_mask_overlapping_partition_forward_;
   std::vector<receive_array> receive_mask_overlapping_partition_backward_;
};

template<typename FMC>
inline void LP<FMC>::construct_factor_partition()
{
    if(factor_partition_valid_) { return; }
    factor_partition_valid_ = true;

    SortFactors();

    union_find uf(f_.size());
    for(auto p : partition_graph) {
        const auto i = factor_address_to_index_[p[0]];
        const auto j = factor_address_to_index_[p[1]];
        uf.merge(i,j);
    }
    auto contiguous_ids = uf.get_contiguous_ids();
    std::vector<INDEX> partition_size(uf.count(),0);
    for(std::size_t i=0; i<contiguous_ids.size(); ++i) {
        const std::size_t id = contiguous_ids[ uf.find(i) ];
        if(f_[i]->FactorUpdated()) {
            partition_size[id]++;
        } 
    }

    // filter out partitions with size zero
    std::vector<INDEX> partition_size_filtered;
    std::vector<INDEX> contiguous_id_to_partition_id(contiguous_ids.size(),std::numeric_limits<INDEX>::max());
    for(INDEX i=0; i<partition_size.size(); ++i) {
        if(partition_size[i]>0) {
            contiguous_id_to_partition_id[i] = partition_size_filtered.size();
            partition_size_filtered.push_back( partition_size[i] );
        }
    }
    //partition_size.erase(std::remove(partition_size.begin(), partition_size.end(), 0), partition_size.end()); 

    factor_partition_ = two_dim_variable_array<FactorTypeAdapter*>(partition_size_filtered);

    // populate factor_partition
    std::fill(partition_size_filtered.begin(), partition_size_filtered.end(), 0);
    for(std::size_t i=0; i<f_.size(); ++i) {
        const std::size_t id = contiguous_id_to_partition_id[ contiguous_ids[ uf.find(i) ] ];
        if(f_[i]->FactorUpdated()) {
            factor_partition_[id][ partition_size_filtered[id]++ ] = f_[i];
        } 
    }

    // sort factor_partition.
    std::unordered_map<FactorTypeAdapter*, std::size_t> factor_address_to_sorted_position;
    for(std::size_t i=0; i<forwardOrdering_.size(); ++i) {
        factor_address_to_sorted_position.insert( {forwardOrdering_[i], i} );
    }
    for(std::size_t i=0; i<factor_partition_.size(); ++i) {
        std::vector<std::pair<std::size_t,FactorTypeAdapter*>> sorted_indices; // sorted index, number in partition
        sorted_indices.reserve(factor_partition_[i].size());
        for(std::size_t j=0; j<factor_partition_[i].size(); ++j) {
            auto* f = factor_partition_[i][j];
            assert(factor_address_to_sorted_position.find(f) != factor_address_to_sorted_position.end());
            const std::size_t idx = factor_address_to_sorted_position[f];
            sorted_indices.push_back( {idx, f} );
        }
        std::sort(sorted_indices.begin(), sorted_indices.end(), [](const auto a, const auto b) { return std::get<0>(a) < std::get<0>(a); });
        for(std::size_t j=0; j<factor_partition_[i].size(); ++j) {
            factor_partition_[i][j] = std::get<1>(sorted_indices[j]);
        } 
    }

    // compute weights and receive masks
    omega_partition_forward_.resize(factor_partition_.size());
    omega_partition_backward_.resize(factor_partition_.size());
    receive_mask_partition_forward_.resize(factor_partition_.size());
    receive_mask_partition_backward_.resize(factor_partition_.size());
    for(std::size_t i=0; i<factor_partition_.size(); ++i) {
        std::tie(omega_partition_forward_[i], receive_mask_partition_forward_[i]) = compute_anisotropic_weights( factor_partition_[i].begin(), factor_partition_[i].end(), 0.0);
        std::tie(omega_partition_backward_[i], receive_mask_partition_backward_[i]) = compute_anisotropic_weights( factor_partition_[i].rbegin(), factor_partition_[i].rend(), 0.0);
    } 

    construct_forward_pushing_weights(factor_partition_.begin(), factor_partition_.end(), omega_partition_forward_pass_push_forward_, receive_mask_partition_forward_pass_push_forward_, 
            [](const auto& partition) { return std::make_tuple(partition.begin(), partition.end()); }
            );

    construct_forward_pushing_weights(factor_partition_.begin(), factor_partition_.end(), omega_partition_backward_pass_push_forward_, receive_mask_partition_backward_pass_push_forward_,
            [](const auto& partition) { return std::make_tuple(partition.rbegin(), partition.rend()); }
            );

    construct_forward_pushing_weights(factor_partition_.rbegin(), factor_partition_.rend(), omega_partition_forward_pass_push_backward_, receive_mask_partition_forward_pass_push_backward_, 
            [](const auto& partition) { return std::make_tuple(partition.begin(), partition.end()); }
            );

    construct_forward_pushing_weights(factor_partition_.rbegin(), factor_partition_.rend(), omega_partition_backward_pass_push_backward_, receive_mask_partition_backward_pass_push_backward_, 
            [](const auto& partition) { return std::make_tuple(partition.rbegin(), partition.rend()); }
            );


    omega_partition_forward_pass_push_.resize(factor_partition_.size()-1);
    receive_mask_partition_forward_pass_push_.resize(factor_partition_.size()-1);
    for(std::size_t i=0; i<factor_partition_.size()-1; ++i) {
        auto f = concatenate_factors(factor_partition_[i].begin(), factor_partition_[i].end(), factor_partition_[i+1].rbegin(), factor_partition_[i+1].rend());
        ComputeAnisotropicWeights( f.begin(), f.end(), omega_partition_forward_pass_push_[i], receive_mask_partition_forward_pass_push_[i], 0.0); 
    }

    omega_partition_backward_pass_push_.resize(factor_partition_.size()-1);
    receive_mask_partition_backward_pass_push_.resize(factor_partition_.size()-1);
    for(std::size_t ri=0; ri<factor_partition_.size()-1; ++ri) {
        const std::size_t i = factor_partition_.size() - ri - 1;
        auto f = concatenate_factors(factor_partition_[i].begin(), factor_partition_[i].end(), factor_partition_[i-1].rbegin(), factor_partition_[i-1].rend());
        ComputeAnisotropicWeights( f.begin(), f.end(), omega_partition_backward_pass_push_[ri], receive_mask_partition_backward_pass_push_[ri], 0.0); 
    }
}

template<typename FMC>
inline void LP<FMC>::construct_overlapping_factor_partition()
{
    construct_factor_partition();
    if(overlapping_factor_partition_valid_) { return; }
    overlapping_factor_partition_valid_ = true;

    omega_overlapping_partition_forward_.resize(factor_partition_.size()-1);
    omega_overlapping_partition_backward_.resize(factor_partition_.size()-1);
    receive_mask_overlapping_partition_forward_.resize(factor_partition_.size()-1);
    receive_mask_overlapping_partition_backward_.resize(factor_partition_.size()-1);
    for(std::size_t i=0; i<factor_partition_.size()-1; ++i) {
        auto f_forward = concatenate_factors(factor_partition_[i].begin(), factor_partition_[i].end(), factor_partition_[i+1].rbegin(), factor_partition_[i+1].rend());
        ComputeAnisotropicWeights( f_forward.begin(), f_forward.end(), omega_overlapping_partition_forward_[i], receive_mask_overlapping_partition_forward_[i], 0.0);

        auto f_backward = concatenate_factors(factor_partition_[i+1].begin(), factor_partition_[i+1].end(), factor_partition_[i].rbegin(), factor_partition_[i].rend());
        ComputeAnisotropicWeights( f_backward.begin(), f_backward.end(), omega_overlapping_partition_backward_[i], receive_mask_overlapping_partition_backward_[i], 0.0);
    } 
}

template<typename FMC>
template<typename PARTITION_ITERATOR, typename INTRA_PARTITION_FACTOR_ITERATOR>
inline void LP<FMC>::construct_forward_pushing_weights(PARTITION_ITERATOR partition_begin, PARTITION_ITERATOR partition_end, std::vector<weight_array>& omega_partition, std::vector<receive_array>& receive_mask_partition, INTRA_PARTITION_FACTOR_ITERATOR factor_iterator_getter)
{
    std::unordered_map<FactorTypeAdapter*, std::size_t> factor_address_to_partition;
    factor_address_to_partition.reserve(f_.size());

    for(auto partition_it=partition_begin; partition_it!=partition_end; ++partition_it) {

        auto [factor_begin, factor_end] = factor_iterator_getter(*partition_it);

        const auto partition_number = std::distance(partition_begin, partition_it);
        for(auto factor_it=factor_begin; factor_it!=factor_end; ++factor_it) {
            assert(factor_address_to_partition.count(*factor_it) == 0);
            factor_address_to_partition.insert( {*factor_it, partition_number} );
        }
    }

    const auto no_partitions = std::distance(partition_begin, partition_end);
    omega_partition.reserve(no_partitions);
    receive_mask_partition.reserve(no_partitions);

    for(auto partition_it=partition_begin; partition_it!=partition_end; ++partition_it) {

        auto [factor_begin, factor_end] = factor_iterator_getter(*partition_it);

        std::unordered_map<FactorTypeAdapter*, std::size_t> factor_address_to_partition_position;
        for(auto factor_it=factor_begin; factor_it!=factor_end; ++factor_it) {
            factor_address_to_partition_position.insert( {*factor_it, factor_address_to_partition_position.size()} );
        }

        const auto partition_number = std::distance(partition_begin, partition_it);
        omega_partition.push_back(allocate_omega(factor_begin, factor_end));
        auto& omega = omega_partition.back();
        receive_mask_partition.push_back(allocate_receive_mask(factor_begin, factor_end));
        auto& receive_mask = receive_mask_partition.back();

        std::size_t c = 0;
        for(auto factor_it=factor_begin; factor_it!=factor_end; ++factor_it) {
            if((*factor_it)->FactorUpdated()) {
                std::size_t k_send = 0;
                std::size_t k_receive = 0;
                for(auto m : (*factor_it)->get_messages()) {
                    if(m.sends_to_adjacent_factor) {
                        if(factor_address_to_partition.count(m.adjacent_factor)) {
                            const auto adjacent_factor_partition = factor_address_to_partition[m.adjacent_factor];
                            if(adjacent_factor_partition >= partition_number) {
                                omega[c][k_send] = 1.0;
                            } else {
                                omega[c][k_send] = 0.0;
                            }
                        } else {
                        //    if(factor_address_to_partition_position.count(m.adjacent_factor) && factor_address_to_partition_position[m.adjacent_factor] > std::distance(factor_begin, factor_it) ) {
                        //        omega[c][k_send] = 1.0;
                        //    } else {
                                omega[c][k_send] = 0.0; 
                        //    }
                        }
                        ++k_send;
                    }
                    if(m.receives_from_adjacent_factor) {
                        const auto adjacent_factor_partition = factor_address_to_partition[m.adjacent_factor];
                        if(adjacent_factor_partition <= partition_number) {
                            receive_mask[c][k_receive] = 1.0;
                        } else {
                            receive_mask[c][k_receive] = 0.0;
                        }
                        ++k_receive;
                    }
                } 
                assert(k_receive == (*factor_it)->no_receive_messages());
                assert(k_send == (*factor_it)->no_send_messages());

                // reweight omega
                const auto omega_sum = std::accumulate(omega[c].begin(), omega[c].end(), 0.0);
                if(omega_sum > 0) {
                    for(auto& x : omega[c]) {
                        x *= 1.0/omega_sum;
                    }
                }
                assert(std::accumulate(omega[c].begin(), omega[c].end(), 0.0) <= 1.0 + eps);
                ++c;
            }
        }
    } 
}

template<typename FMC> 
void LP<FMC>::compute_partition_pass(const std::size_t no_passes)
{
    construct_factor_partition();
    for(std::size_t i=0; i<factor_partition_.size(); ++i) {
        for(std::size_t iter=0; iter<no_passes; ++iter) {
            ComputePass(factor_partition_[i].begin(), factor_partition_[i].end(), omega_partition_forward_[i].begin(), receive_mask_partition_forward_[i].begin());
            ComputePass(factor_partition_[i].rbegin(), factor_partition_[i].rend(), omega_partition_backward_[i].begin(), receive_mask_partition_backward_[i].begin());
        } 
        // push all messages forward
        if(i < factor_partition_.size()-1) {
            auto f = concatenate_factors(factor_partition_[i].begin(), factor_partition_[i].end(), factor_partition_[i+1].rbegin(), factor_partition_[i+1].rend());
            ComputePass(f.begin(), f.end(), omega_partition_forward_pass_push_[i].begin(), receive_mask_partition_forward_pass_push_[i].begin());
        }
        //ComputePass(factor_partition_[i].begin(), factor_partition_[i].end(), omega_partition_forward_pass_push_forward_[i].begin(), receive_mask_partition_forward_pass_push_forward_[i].begin());
        //ComputePass(factor_partition_[i].rbegin(), factor_partition_[i].rend(), omega_partition_backward_pass_push_forward_[i].begin(), receive_mask_partition_backward_pass_push_forward_[i].begin());
    } 

    for(std::size_t ri=0; ri<factor_partition_.size(); ++ri) {
        const std::size_t i = factor_partition_.size() - ri - 1;
        for(std::size_t iter=0; iter<no_passes; ++iter) {
            ComputePass(factor_partition_[i].begin(), factor_partition_[i].end(), omega_partition_forward_[i].begin(), receive_mask_partition_forward_[i].begin());
            ComputePass(factor_partition_[i].rbegin(), factor_partition_[i].rend(), omega_partition_backward_[i].begin(), receive_mask_partition_backward_[i].begin());
        } 
        // push all messages backward
        if(i != 0) {
            auto f = concatenate_factors(factor_partition_[i].begin(), factor_partition_[i].end(), factor_partition_[i-1].rbegin(), factor_partition_[i-1].rend());
            ComputePass(f.begin(), f.end(), omega_partition_backward_pass_push_[ri].begin(), receive_mask_partition_backward_pass_push_[ri].begin());
        }
        //ComputePass(factor_partition_[i].begin(), factor_partition_[i].end(), omega_partition_forward_pass_push_backward_[ri].begin(), receive_mask_partition_forward_pass_push_backward_[ri].begin());
        //ComputePass(factor_partition_[i].rbegin(), factor_partition_[i].rend(), omega_partition_backward_pass_push_backward_[ri].begin(), receive_mask_partition_backward_pass_push_backward_[ri].begin());
    } 
}

template<typename FMC> 
void LP<FMC>::compute_overlapping_partition_pass(const std::size_t no_passes)
{
    construct_overlapping_factor_partition();

/*
    auto forward_pass = [&](const std::size_t begin, const std::size_t end) 
    {
        for(std::size_t i=begin; i<end; ++i) {
            auto f_forward = concatenate_factors(factor_partition_[i].begin(), factor_partition_[i].end(), factor_partition_[i+1].rbegin(), factor_partition_[i+1].rend());
            auto f_backward = concatenate_factors(factor_partition_[i+1].begin(), factor_partition_[i+1].end(), factor_partition_[i].rbegin(), factor_partition_[i].rend());
            for(std::size_t iter=0; iter<no_passes; ++iter) {
                ComputePass( f_forward.begin(), f_forward.end(), omega_overlapping_partition_forward_[i].begin(), receive_mask_overlapping_partition_forward_[i].begin());
                ComputePass( f_backward.begin(), f_backward.end(), omega_overlapping_partition_backward_[i].begin(), receive_mask_overlapping_partition_backward_[i].begin());
            }
            ComputePass( f_forward.begin(), f_forward.end(), omega_overlapping_partition_forward_[i].begin(), receive_mask_overlapping_partition_forward_[i].begin());
        }
    };

    auto backward_pass = [&](const std::size_t begin, const std::size_t end)
    {
        for(std::size_t ri=begin+1; ri<end+1; ++ri) {
            const std::size_t i = factor_partition_.size() - ri - 1;
            auto f_forward = concatenate_factors(factor_partition_[i].begin(), factor_partition_[i].end(), factor_partition_[i+1].rbegin(), factor_partition_[i+1].rend());
            auto f_backward = concatenate_factors(factor_partition_[i+1].begin(), factor_partition_[i+1].end(), factor_partition_[i].rbegin(), factor_partition_[i].rend());

            for(std::size_t iter=0; iter<no_passes; ++iter) {
                ComputePass( f_backward.begin(), f_backward.end(), omega_overlapping_partition_backward_[i].begin(), receive_mask_overlapping_partition_backward_[i].begin());
                ComputePass( f_forward.begin(), f_forward.end(), omega_overlapping_partition_forward_[i].begin(), receive_mask_overlapping_partition_forward_[i].begin());
            }
            ComputePass( f_backward.begin(), f_backward.end(), omega_overlapping_partition_backward_[i].begin(), receive_mask_overlapping_partition_backward_[i].begin());
        } 
    };

    const auto no_threads = 1;
    const auto loop_size = factor_partition_.size()-1;
    std::vector<std::thread> threads(no_threads);
    const int grainsize = loop_size / no_threads;

    auto work_iter = 0;
    for(auto it = std::begin(threads); it != std::end(threads) - 1; ++it) {
        *it = std::thread(forward_pass, work_iter, work_iter + grainsize);
        work_iter += grainsize;
    }
    threads.back() = std::thread(forward_pass, work_iter, factor_partition_.size()-1);

    for(auto&& i : threads) {
        i.join();
    }

    work_iter = 0;
    for(auto it = std::begin(threads); it != std::end(threads) - 1; ++it) {
        *it = std::thread(backward_pass, work_iter, work_iter + grainsize);
        work_iter += grainsize;
    }
    threads.back() = std::thread(backward_pass, work_iter, factor_partition_.size()-1);

    for(auto&& i : threads) {
        i.join();
    }


    return;
*/

    for(std::size_t i=0; i<factor_partition_.size()-1; ++i) {
        auto f_forward = concatenate_factors(factor_partition_[i].begin(), factor_partition_[i].end(), factor_partition_[i+1].rbegin(), factor_partition_[i+1].rend());
        auto f_backward = concatenate_factors(factor_partition_[i+1].begin(), factor_partition_[i+1].end(), factor_partition_[i].rbegin(), factor_partition_[i].rend());
        for(std::size_t iter=0; iter<no_passes; ++iter) {
            ComputePass( f_forward.begin(), f_forward.end(), omega_overlapping_partition_forward_[i].begin(), receive_mask_overlapping_partition_forward_[i].begin());
            ComputePass( f_backward.begin(), f_backward.end(), omega_overlapping_partition_backward_[i].begin(), receive_mask_overlapping_partition_backward_[i].begin());
        }
        ComputePass( f_forward.begin(), f_forward.end(), omega_overlapping_partition_forward_[i].begin(), receive_mask_overlapping_partition_forward_[i].begin());
    } 

    for(std::size_t ri=1; ri<factor_partition_.size(); ++ri) {
        const std::size_t i = factor_partition_.size() - ri - 1;
        auto f_forward = concatenate_factors(factor_partition_[i].begin(), factor_partition_[i].end(), factor_partition_[i+1].rbegin(), factor_partition_[i+1].rend());
        auto f_backward = concatenate_factors(factor_partition_[i+1].begin(), factor_partition_[i+1].end(), factor_partition_[i].rbegin(), factor_partition_[i].rend());

        for(std::size_t iter=0; iter<no_passes; ++iter) {
            ComputePass( f_backward.begin(), f_backward.end(), omega_overlapping_partition_backward_[i].begin(), receive_mask_overlapping_partition_backward_[i].begin());
            ComputePass( f_forward.begin(), f_forward.end(), omega_overlapping_partition_forward_[i].begin(), receive_mask_overlapping_partition_forward_[i].begin());
        }
        ComputePass( f_backward.begin(), f_backward.end(), omega_overlapping_partition_backward_[i].begin(), receive_mask_overlapping_partition_backward_[i].begin());
    }
}

template<typename FMC>
template<typename ITERATOR_1, typename ITERATOR_2>
std::vector<FactorTypeAdapter*> LP<FMC>::concatenate_factors(ITERATOR_1 f1_begin, ITERATOR_1 f1_end, ITERATOR_2 f2_begin, ITERATOR_2 f2_end)
{
    std::vector<FactorTypeAdapter*> f;
    f.reserve(std::distance(f1_begin, f1_end) + std::distance(f2_begin, f2_end));
    std::copy(f1_begin, f1_end, std::back_inserter(f));
    std::copy(f2_begin, f2_end, std::back_inserter(f)); 
    return f; 
}

} // namespace LPMP

#endif // LPMP_FACTOR_PARTITIONS_STORAGE_HXX

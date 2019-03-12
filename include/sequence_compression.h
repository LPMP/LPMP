#pragma once

#include <vector>
#include <limits>
#include <cassert>

namespace LPMP {

    class sequence_compression {
        public:
            constexpr static std::size_t index_not_present = std::numeric_limits<std::size_t>::max();

            sequence_compression(const std::size_t N)
                : orig_to_compressed_index_(N, index_not_present)
            {}

            void reset()
            {
                for(const std::size_t orig_index : compressed_to_orig_index_)
                    orig_to_compressed_index_[orig_index] = index_not_present;
                compressed_to_orig_index_.clear();
                for(const auto v : orig_to_compressed_index_) {
                    assert(v == index_not_present);
                } 
            }
            std::size_t add_index(const std::size_t index)
            {
                assert(index < orig_to_compressed_index_.size());
                if(orig_to_compressed_index_[index] == index_not_present) {
                    orig_to_compressed_index_[index] = compressed_to_orig_index_.size();
                    compressed_to_orig_index_.push_back(index);
                }
                return compressed_to_orig_index_.size()-1;
            }

    std::size_t no_compressed_indices() const { return compressed_to_orig_index_.size(); }
            std::size_t no_orig_indices() const { return orig_to_compressed_index_.size(); }
            std::size_t orig_to_compressed_index(const std::size_t orig_index) const
            {
                assert(orig_index < orig_to_compressed_index_.size());
                assert(orig_to_compressed_index_[orig_index] != index_not_present);
                assert(orig_to_compressed_index_[orig_index] < orig_to_compressed_index_.size());
                return orig_to_compressed_index_[orig_index];
            }

            std::size_t compressed_to_orig_index(const std::size_t compressed_index) const
            {
                assert(compressed_index < compressed_to_orig_index_.size());
                assert(compressed_index == orig_to_compressed_index_[ compressed_to_orig_index_[compressed_index] ]);
                return compressed_to_orig_index_[compressed_index]; 
            }

        private:
            std::vector<std::size_t> orig_to_compressed_index_;
            std::vector<std::size_t> compressed_to_orig_index_; 
    };
} // namespace LPMP

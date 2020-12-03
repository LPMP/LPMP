#pragma once

#include <vector>

namespace LPMP {

    class permutation : public std::vector<std::size_t> {
        public:
            using std::vector<std::size_t>::vector;
            bool is_permutation() const;
            permutation inverse_permutation() const;
            template<typename ITERATOR>
                std::vector<typename std::iterator_traits<ITERATOR>::value_type> permute(ITERATOR begin, ITERATOR end) const;
    };

    template<typename ITERATOR>
        bool is_permutation(ITERATOR begin, ITERATOR end)
        {
            std::vector<char> nr_taken(std::distance(begin, end), 0);
            for(auto it=begin; it!=end; ++it) {
                if(*it >= std::distance(begin,end))
                    return false;
                if(nr_taken[*it] != 0)
                    return false;
                nr_taken[*it] = 1; 
            }
            return true;
        }

    template<typename ITERATOR>
        permutation inverse_permutation(ITERATOR begin, ITERATOR end)
        {
            permutation inverse_perm(std::distance(begin, end));
            for(std::size_t i=0; i<std::distance(begin, end); ++i) 
                inverse_perm[(*begin+i)] = i;
            return inverse_perm;
        }

    inline bool permutation::is_permutation() const
    {
        return LPMP::is_permutation(this->begin(), this->end());
    }

    inline permutation permutation::inverse_permutation() const
    {
        assert(is_permutation());
        return LPMP::inverse_permutation(this->begin(), this->end());
    }

    template<typename ITERATOR>
        std::vector<typename std::iterator_traits<ITERATOR>::value_type> permutation::permute(ITERATOR begin, ITERATOR end) const
        {
            assert(is_permutation());
            assert(std::distance(begin,end) == this->size());
            std::vector<typename std::iterator_traits<ITERATOR>::value_type> result;
            result.reserve(this->size());
            for(std::size_t i=0; i<this->size(); ++i) {
                result.push_back( *(begin+(*this)[i]) );
            }

            return result;
        }
}

#pragma once

#include <initializer_list>
#include "memory_allocator.hxx"
//#include "serialization.hxx"
#include "config.hxx"
#include "help_functions.hxx"
//#include "cereal/archives/binary.hpp"
//#include <heaplayers>
//#include <heaps/special/obstackheap.h>
//#include "aligned_allocator.h"

namespace LPMP {


    // fixed size vector allocated from block allocator
    // possibly holding size explicitly is not needed: It is held by allocator as well

    // primitive expression templates for vector
    template<typename T, typename E>
        class vector_expression {
            public:
                T operator[](const std::size_t i) const { return static_cast<E const&>(*this)[i]; }
                T operator()(const std::size_t i1, const std::size_t i2) const { return static_cast<E const&>(*this)(i1,i2); }
                T operator()(const std::size_t i1, const std::size_t i2, const std::size_t i3) const { return static_cast<E const&>(*this)(i1,i2,i3); }
                std::size_t size() const { return static_cast<E const&>(*this).size(); }
                std::size_t dim1() const { return static_cast<E const&>(*this).dim1(); }
                std::size_t dim2() const { return static_cast<E const&>(*this).dim2(); }
                std::size_t dim3() const { return static_cast<E const&>(*this).dim3(); }

                E& operator()() { return static_cast<E&>(*this); }
                const E& operator()() const { return static_cast<const E&>(*this); }

                friend std::ostream& operator<<(std::ostream& os, const vector_expression<T,E>& v) {
                    for(std::size_t i=0; i<v.size(); ++i) {
                        os << v[i] << " ";
                    }
                    os << "\n";
                    return os;
                }
        };

    template<typename T, typename E>
        class matrix_expression {
            public:
                T operator[](const std::size_t i) const { return static_cast<E const&>(*this)[i]; }
                std::size_t size() const { return static_cast<E const&>(*this).size(); }

                T operator()(const std::size_t i1, const std::size_t i2) const { return static_cast<E const&>(*this)(i1,i2); }
                std::size_t dim1() const { return static_cast<E const&>(*this).dim1(); }
                std::size_t dim2() const { return static_cast<E const&>(*this).dim2(); }

                E& operator()() { return static_cast<E&>(*this); }
                const E& operator()() const { return static_cast<const E&>(*this); }

                friend std::ostream& operator<<(std::ostream& os, const matrix_expression<T,E>& m) {
                    for(std::size_t i=0; i<m.dim1(); ++i) {
                        for(std::size_t j=0; j<m.dim2(); ++j) {
                            os << m(i,j) << " ";
                        }
                        os << "\n";
                    }
                    return os;
                }
        };

    template<typename T, typename E>
        class tensor3_expression {
            public:
                T operator[](const std::size_t i) const { return static_cast<E const&>(*this)[i]; }
                std::size_t size() const { return static_cast<E const&>(*this).size(); }

                T operator()(const std::size_t i1, const std::size_t i2) const { return static_cast<E const&>(*this)(i1,i2); }
                T operator()(const std::size_t i1, const std::size_t i2, const std::size_t i3) const { return static_cast<E const&>(*this)(i1,i2,i3); }
                std::size_t dim1() const { return static_cast<E const&>(*this).dim1(); }
                std::size_t dim2() const { return static_cast<E const&>(*this).dim2(); }
                std::size_t dim3() const { return static_cast<E const&>(*this).dim3(); }

                E& operator()() { return static_cast<E&>(*this); }
                const E& operator()() const { return static_cast<const E&>(*this); }
        };


    // possibly also support different allocators: a pure stack allocator without block might be a good choice as well for short-lived memory
    template<typename T=REAL>
        class vector : public vector_expression<T,vector<T>> {
            public:

                void fill()
                {
                    std::fill(begin_, end_, T{});
                }
                template<typename ITERATOR>
                    vector(ITERATOR begin, ITERATOR end)
                    {
                        const std::size_t size = std::distance(begin,end);
                        assert(size > 0);
                        //begin_ = (T*) std::aligned_alloc(32, size*sizeof(T));
                        begin_ = new T[size];
                        assert(begin_ != nullptr);
                        end_ = begin_ + size;
                        for(auto it=this->begin(); begin!=end; ++begin, ++it) {
                            (*it) = *begin;
                        }
                    }

                vector(const std::size_t size) 
                {
                    //begin_ = (T*) std::aligned_alloc(32, size*sizeof(T));
                    begin_ = new T[size];
                    assert(size > 0);
                    assert(begin_ != nullptr);
                    end_ = begin_ + size;
                    fill();
                }
                vector(const std::size_t size, const T value) 
                    : vector(size)
                {
                    assert(size > 0);
                    std::fill(begin_, end_, value);
                }
                vector()
                    : begin_(nullptr),
                    end_(nullptr)
            {}
                ~vector() {
                    if(begin_ != nullptr) {
                        //global_real_block_arena_array[stack_allocator_index].deallocate((void*)begin_);
                        std::free(begin_);
                        //allocator::get().free(begin_);
                    }
                    static_assert(sizeof(T) % sizeof(int) == 0,"");
                }
                vector(const vector& o)  
                    : vector(o.size())
                {
                    assert(begin_ != nullptr);
                    auto it = begin_;
                    for(auto o_it = o.begin(); o_it!=o.end(); ++it, ++o_it) { *it = *o_it; }
                }

                vector(vector&& o)
                {
                    assert(begin_ != o.begin_ && end_ != o.end_);
                    std::swap(begin_, o.begin_);
                    std::swap(end_, o.end_);
                }

                vector& operator=(const vector<T>& o)
                {
                    if(size() != o.size()) {
                        vector copy(o.size());
                        std::swap(begin_, copy.begin_);
                        std::swap(end_, copy.begin_);
                    }
                    assert(size() == o.size());
                    for(std::size_t i=0; i<o.size(); ++i) { 
                        (*this)[i] = o[i]; 
                    } 
                    return *this;
                }

                vector& operator=(vector&& o)
                {
                    std::swap(begin_, o.begin_);
                    std::swap(end_, o.end_);
                    return *this;
                }

                vector& operator+=(const vector<T>& o)
                {
                    static_assert(std::is_same<T,double>::value,"");
                    assert(size() == o.size());
                    for(std::size_t i=0; i<size(); i++) {
                        (*this)[i] += o[i];
                    }
                    return *this;
                }

                void share(const vector& o)
                {
                    // memory leak here!
                    assert(false);
                    begin_ = o.begin_;
                    end_ = o.end_; 
                }

                template<typename E>
                    bool operator==(const vector_expression<T,E>& o) {
                        assert(size() == o.size());
                        for(std::size_t i=0; i<o.size(); ++i) { 
                            if((*this)[i] != o[i]) { 
                                return false;
                            }
                        }
                        return true;
                    }
                template<typename E>
                    void operator=(const vector_expression<T,E>& o) {
                        assert(size() == o.size());
                        for(std::size_t i=0; i<o.size(); ++i) { 
                            (*this)[i] = o[i]; 
                        }
                    }
                template<typename E>
                    void operator-=(const vector_expression<T,E>& o) {
                        assert(size() == o.size());
                        for(std::size_t i=0; i<o.size(); ++i) { 
                            (*this)[i] -= o[i]; } 
                    }
                template<typename E>
                    void operator+=(const vector_expression<T,E>& o) {
                        assert(size() == o.size());
                        for(std::size_t i=0; i<o.size(); ++i) { 
                            (*this)[i] += o[i]; } 
                    }

                void operator+=(const T x) {
                    for(std::size_t i=0; i<this->size(); ++i) { 
                        (*this)[i] += x;
                    } 
                }

                void operator-=(const T x) {
                    for(std::size_t i=0; i<this->size(); ++i) { 
                        (*this)[i] -= x;
                    } 
                }


                // force construction from expression template
                template<typename E>
                    vector(vector_expression<T,E>& v) : vector(v.size())
                {
                    for(std::size_t i=0; v.size(); ++i) {
                        (*this)[i] = v[i];
                    }
                }

                std::size_t size() const { return end_ - begin_; }

                const T& operator[](const std::size_t i) const {
                    assert(i<size());
                    //assert(!std::isnan(begin_[i]));
                    return begin_[i];
                }
                T& operator[](const std::size_t i) {
                    assert(i<size());
                    //assert(!std::isnan(begin_[i]));
                    return begin_[i];
                }
                using iterator = T*;
                T* begin() const { return begin_; }
                T* end() const { return end_; }

                T back() const { assert(size() > 0); return *(end_-1); }
                T& back() { assert(size() > 0); return *(end_-1); }

                template<typename ARCHIVE>
                    void serialize(ARCHIVE& ar)
                    {
                        assert(false);
                        //ar( cereal::binary_data( begin_, sizeof(T)*size()) );
                        //ar( binary_data<REAL>( begin_, sizeof(T)*size()) );
                    }

                void prefetch() const { }

                // minimum operation with simd instructions (when T is float, double or integer)
                T min() const
                {
                    //assert((std::size_t(begin_) % 32) == 0);
                    return *std::min_element(begin(), end());
                }

                T min_except(const std::size_t i) const
                {
                    assert(i < this->size());
                    const auto val = (*this)[i];
                    begin_[i] = std::numeric_limits<REAL>::infinity();
                    const auto min_val = this->min();
                    begin_[i] = val;
                    return min_val;
                }

                std::array<T,2> two_min() const
                {
                    assert(size() >= 2);
                    REAL min_val = std::numeric_limits<REAL>::infinity();
                    REAL second_min_val = std::numeric_limits<REAL>::infinity();
                    for(size_t i=0; i<this->size(); ++i)
                    {
                        const REAL tmp_min = std::min((*this)[i], min_val);
                        const REAL tmp_max = std::max((*this)[i], min_val);
                        min_val = tmp_min;
                        second_min_val = std::min(tmp_max, second_min_val);
                    }

                    return {min_val, second_min_val};
                }

    private:
    T* begin_ = nullptr;
    T* end_ = nullptr;
};

// vector that caches minimum value
namespace min_vector_update_handle {
    template<typename T>
        class type {
            public:
                type(T& val, T& min_value, bool& min_value_valid)
                    : val_(val),
                    min_value_(min_value),
                    min_value_valid_(min_value_valid)
            {}

                template<typename T2>
                    friend void swap(type<T2>& t1, type<T2>& t2);

                operator T&() const { return val_; }

                type& operator=(const T& x) { set_value(x); return *this; }
                type& operator=(const type& o) { set_value(o.val_); return *this; }
                type& operator*(const T& x) { set_value(val_*x); return *this; }
                type& operator/(const T& x) { set_value(val_/x); return *this; } 
                type& operator-=(const T& x) { change_value(-x); return *this; }
                type& operator+=(const T& x) { change_value(x); return *this; }

            private:
                void update_value(const T& x)
                {
                    if(x <= 0.0) {
                        val_ += x;
                        min_value_ = std::min(min_value_, val_); 
                    } else {
                        if(val_ == min_value_) {
                            min_value_valid_ = false; 
                        }
                        val_ += x;
                    }
                }

                void set_value(const T& x)
                {
                    if(x < min_value_) {
                        min_value_ = x;
                    } else if(val_ == min_value_) {
                        min_value_valid_ = false; 
                    }
                    val_ = x;
                }

                T& val_;
                T& min_value_;
                bool& min_value_valid_;
        }; 

    template<typename T>
        void swap(type<T>& t1, type<T>& t2) 
        { 
            std::swap(t1.ptr_, t2.ptr_); 
        } 
}

template<typename T=REAL>
class min_vector : public vector_expression<T,min_vector<T>> {
    public:
        // TODO: perfect forwarding
        template<typename... ARGS>
            min_vector(ARGS... args)
            : vec_(args...)
            {}

        auto size() const { return vec_.size(); }

        template<typename E>
            void operator=(const vector_expression<T,E>& o) {
                assert(size() == o.size());
                min_value_ = o[0];
                for(std::size_t i=0; i<o.size(); ++i) { 
                    (*this)[i] = o[i]; 
                    min_value_ = std::min(min_value_, o[i]);
                }
                min_value_valid_ = true;
            }
        template<typename E>
            void operator-=(const vector_expression<T,E>& o) {
                assert(size() == o.size());
                min_value_ = o[0];
                for(std::size_t i=0; i<o.size(); ++i) { 
                    (*this)[i] -= o[i]; 
                    min_value_ = std::min(min_value_, (*this)[i]);
                } 
                min_value_valid_ = true;
            }
        template<typename E>
            void operator+=(const vector_expression<T,E>& o) {
                assert(size() == o.size());
                min_value_ = o[0];
                for(std::size_t i=0; i<o.size(); ++i) { 
                    (*this)[i] += o[i]; 
                    min_value_ = std::min(min_value_, (*this)[i]);
                } 
                min_value_valid_ = true;
            }

        T min() const
        {
            if(!min_value_valid_) {
                min_value_ = vec_.min();
                min_value_valid_ = true;
            }

            assert(min_value_ == vec_.min());
            return min_value_;
        }

        T min_except(const std::size_t i) const
        {
            if(min_value_valid_ && vec_[i] > min_value_) {
                return min_value_;
            } else {
                const T v = vec_.min_except(i);
                min_value_ = std::min(v, vec_[i]);
                min_value_valid_ = true;
                return v; 
            }
        }

        std::array<T,2> two_min() const
        {
            return vec_.two_min();
        }

        T& operator[](const std::size_t i) {
            assert(i < this->size());
            return min_vector_update_handle::type<T>( vec_[i], min_value_, min_value_valid_ );
        }
        const T& operator[](const std::size_t i) const { return vec_[i]; }

        // does not work properly: operator* expects reference
        /*
           class iterator {
           public:
           using difference_type = std::ptrdiff_t;
           using value_type  = min_vector_update_handle::type<T>;
           using pointer = T*;
           using reference = T&;
           using iterator_category = std::random_access_iterator_tag;

           iterator(T* ptr, T& min_value, bool& min_value_valid)
           : ptr_(ptr),
           min_value_(min_value),
           min_value_valid_(min_value_valid)
           {}

           iterator& operator=(const iterator& o)
           {
           assert(&min_value_ == &o.min_value_);
           ptr_ = o.ptr_; 
           return *this;
           }

           bool operator==(const iterator& o) const { return ptr_ == o.ptr_; }
           bool operator!=(const iterator& o) const { return ptr_ != o.ptr_; }
           iterator& operator++() { ++ptr_; return *this; }
           iterator& operator--() { --ptr_; return *this; }
           iterator operator+(const std::size_t n) const { return iterator(ptr_ + n, min_value_, min_value_valid_); }
           iterator operator-(const std::size_t n) const { return iterator(ptr_ - n, min_value_, min_value_valid_); }
           std::ptrdiff_t operator-(const iterator& o) const { return ptr_ - o.ptr_; }
           std::ptrdiff_t operator+(const iterator& o) const { return ptr_ + o.ptr_; }
           bool operator<(const iterator& o) const { return ptr_ < o.ptr_; }

           min_vector_update_handle::type<T> operator*() { return min_vector_update_handle::type<T>(*ptr_, min_value_, min_value_valid_); }
           min_vector_update_handle::type<T> operator->() { return min_vector_update_handle::type<T>(*ptr_, min_value_, min_value_valid_); }

           private:
           T* ptr_;
           T& min_value_;
           bool& min_value_valid_; 
           };

           iterator begin() { return iterator(vec_.begin(), min_value_, min_value_valid_); }
           iterator end() { return iterator(vec_.end(), min_value_, min_value_valid_); }
           */

        T* begin() { min_value_valid_ = false; return vec_.begin(); }
        T* end() { min_value_valid_ = false; return vec_.end(); }

        min_vector_update_handle::type<T> back() { return min_vector_update_handle::type<T>( vec_.back(), min_value_, min_value_valid_ ); }
        const T& back() const { return vec_.back(); }

        vector<T>& data() { min_value_valid_ = false; return vec_; }

    private:
        vector<T> vec_;
        mutable T min_value_;
        mutable bool min_value_valid_ = false;
};

// additionally store free space so that no allocation is needed for small objects
template<typename T, std::size_t N>
class small_vector : public vector<T>
{
    protected:
        void allocate(const std::size_t s)
        {
            if(s <= N) {
                this->begin_ = local_.begin();
                this->end_ = local_.begin() + s;
                std::fill(this->end(), local_.end(), std::numeric_limits<REAL>::infinity());
            } else {
                vector<T>::allocate(s);
            }
        }

    public:
        template<typename ITERATOR>
            small_vector(ITERATOR begin, ITERATOR end)
            {
                allocate();
            }

        small_vector(const std::size_t size)
        {
            allocate(size);
        }

        small_vector(const std::size_t size, const T val)
        {
            allocate(size);
            std::fill(this->begin(), this->end(), val);
        }

        ~small_vector()
        {
            if(!is_small()) {
                this->deallocate();
            } 
        }

    private:
        bool is_small() const { return this->size() < N; }
        std::array<T,N> local_;
};


template<typename T, std::size_t N>
using array_impl = std::array<T,N>;



template<typename T, std::size_t N>
class array : public vector_expression<T,array<T,N>> {
    public:
        array() {}

        template<typename ITERATOR>
            array(ITERATOR begin, ITERATOR end)
            {
                const std::size_t size = std::distance(begin,end);
                assert(size == N);
                for(auto it=array_.begin(); it!=array_.end(); ++it) {
                    (*it) = *begin;
                }
            }

        array(std::initializer_list<T> l)
            : array(l.begin(), l.end())
        {}

        array(const T value) 
        {
            std::fill(array_.begin(), array_.end(), value);
        }
        array(const array<T,N>& o)  {
            assert(size() == o.size());
            auto it = begin();
            for(auto o_it = o.begin(); o_it!=o.end(); ++it, ++o_it) { *it = *o_it; }
        }
        template<typename E>
            void operator=(const vector_expression<T,E>& o) {
                assert(size() == o.size());
                for(std::size_t i=0; i<o.size(); ++i)
                    (*this)[i] = o[i]; 
            }
        template<typename E>
            void operator-=(const vector_expression<T,E>& o) {
                assert(size() == o.size());
                for(std::size_t i=0; i<o.size(); ++i) { 
                    (*this)[i] -= o[i]; } 
            }
        template<typename E>
            void operator+=(const vector_expression<T,E>& o) {
                assert(size() == o.size());
                for(std::size_t i=0; i<o.size(); ++i) { 
                    (*this)[i] += o[i]; } 
            }


        // force construction from expression template
        template<typename E>
            array(vector_expression<T,E>& v) 
            {
                for(std::size_t i=0; v.size(); ++i) {
                    (*this)[i] = v[i];
                }
            }

        constexpr static std::size_t size() { return N; }

        T operator[](const std::size_t i) const {
            assert(i<size());
            //assert(!std::isnan(array_[i]));
            return array_[i];
        }
        T& operator[](const std::size_t i) {
            assert(i<size());
            return array_[i];
        }
        using iterator = T*;
        auto begin() const { return array_.begin(); }
        auto end() const { return array_.end(); }
        auto begin() { return array_.begin(); }
        auto end() { return array_.end(); }

        template<typename ARCHIVE>
            void serialize(ARCHIVE& ar)
            {
                ar( array_ );
            } 

        T min() const
        {
            static_assert(std::is_same<T,REAL>::value,"");
            return *std::min_element(array_.begin(), array_.end());
        }
    private:
        std::array<T,N> array_; // should only be done for REALs, not generally
};


// matrix is based on vector.
template<typename T=REAL>
class matrix : public matrix_expression<T,matrix<T>> {
    public:
        T& operator[](const std::size_t i) 
        { 
            assert(i<size());
            const std::size_t x1 = i/dim2();
            const std::size_t x2 = i%dim2();
            return (*this)(x1,x2);
        }
        const T operator[](const std::size_t i) const 
        { 
            assert(i<size());
            const std::size_t x1 = i/dim2();
            const std::size_t x2 = i%dim2();
            return (*this)(x1,x2);
        }
        std::size_t size() const { return dim1()*dim2(); }
        template<typename ARCHIVE>
            void serialize(ARCHIVE& ar)
            {
                ar( vec_ );
            }

        class const_iterator {
            public:
                const_iterator(const std::size_t _dim2, const vector<T>& _vec, const std::size_t _i)
                    :
                        i(_i),
                        dim2(_dim2),
                        vec(_vec)
            {}

                void operator++() {
                    ++i; 
                }
                T operator*() const {
                    const std::size_t x1 = i/dim2;
                    const std::size_t x2 = i%dim2;
                    return vec[x1*dim2 + x2]; 
                }
                bool operator==(const const_iterator& o) const {
                    assert(dim2 == o.dim2 && &vec == &o.vec);
                    return i == o.i;
                }
                bool operator!=(const const_iterator& o) const {
                    return !( (*this) == o );
                }
            private:
                std::size_t i;
                const std::size_t dim2;
                const vector<T>& vec;
        };

        class iterator {
            public:
                iterator(const std::size_t _dim2, vector<T>& _vec, const std::size_t _i)
                    :
                        i(_i),
                        dim2(_dim2),
                        vec(_vec)
            {}

                void operator++() {
                    ++i; 
                }
                T& operator*() {
                    const std::size_t x1 = i/dim2;
                    const std::size_t x2 = i%dim2;
                    return vec[x1*dim2 + x2];
                }
                T operator*() const {
                    const std::size_t x1 = i/dim2;
                    const std::size_t x2 = i%dim2;
                    return vec[x1*dim2 + x2]; 
                }
                bool operator==(const iterator& o) const {
                    assert(dim2 == o.dim2 && vec == o.vec);
                    return i == o.i;
                }
                bool operator!=(const iterator& o) const {
                    return !( (*this) == o );
                }
            private:
                std::size_t i;
                const std::size_t dim2;
                vector<T>& vec;
        };

        const_iterator begin() const { return const_iterator(dim2(), dim2(), vec_, 0); }
        const_iterator end() const { return const_iterator(dim2(), dim2(), vec_, size()); }

        iterator begin() { return iterator(dim2(), vec_, 0); }
        iterator end() { return iterator(dim2(), vec_, size()); }

        matrix(const std::size_t d1, const std::size_t d2) : vec_(d1*d2), dim2_(d2) {
            assert(d1 > 0 && d2 > 0);
        }
        matrix(const std::size_t d1, const std::size_t d2, const T val) : vec_(d1*d2), dim2_(d2) {
            assert(d1 > 0 && d2 > 0);
            std::fill(vec_.begin(), vec_.end(), val);
        }
        matrix() : vec_(), dim2_(0) {}
        matrix(const matrix& o) 
            : vec_(o.vec_),
            dim2_(o.dim2_)
    {}
        matrix(matrix&& o) 
            : vec_(std::move(o.vec_)),
            dim2_(o.dim2_)
    {}
        bool operator==(const matrix<T>& o) {
            assert(this->size() == o.size() && o.dim2_ == dim2_);
            return vec_ == o.vec_; 
        }
        void operator=(const matrix<T>& o) {
            if(!(this->vec_.size() == o.vec_.size() && o.dim2_ == dim2_)) {
                matrix copy(o.dim1(), o.dim2());
                std::swap(vec_, copy.vec_);
                std::swap(dim2_, copy.dim2_);
            }
            assert(this->size() == o.size() && o.dim2_ == dim2_);
            vec_ = o.vec_;
        }

        matrix& operator+=(const matrix<T>& o)
        {
            static_assert(std::is_same<T,double>::value,"");
            assert(dim1() == o.dim1());
            assert(dim2() == o.dim2());
            vec_ += o.vec_;
            return *this;
        }

        template<typename E>
            void operator-=(const matrix_expression<T,E>& o) {
                assert(dim1() == o.dim1() && dim2() == o.dim2());
                for(std::size_t i=0; i<o.dim1(); ++i) { 
                    for(std::size_t j=0; j<o.dim2(); ++j) { 
                        (*this)(i,j) -= o(i,j); 
                    } 
                }
            }


        T& operator()(const std::size_t x1, const std::size_t x2) { assert(x1<dim1() && x2<dim2()); return vec_[x1*dim2() + x2]; }
        const T& operator()(const std::size_t x1, const std::size_t x2) const { assert(x1<dim1() && x2<dim2()); return vec_[x1*dim2() + x2]; }
        T& operator()(const std::size_t x1, const std::size_t x2, const std::size_t x3) { assert(x3 == 0); return (*this)(x1,x2); } // sometimes we treat a matrix as a tensor with trivial last dimension
        const T& operator()(const std::size_t x1, const std::size_t x2, const std::size_t x3) const { assert(x3 == 0); return (*this)(x1,x2); } // sometimes we treat a matrix as a tensor with trivial last dimension

        const std::size_t dim1() const { return vec_.size()/dim2(); }
        const std::size_t dim2() const { return dim2_; }

        void transpose() {
            assert(dim1() == dim2());
            for(std::size_t x1=0; x1<dim1(); ++x1) {
                for(std::size_t x2=0; x2<x1; ++x2) {
                    std::swap((*this)(x1,x2), (*this)(x2,x1));
                }
            }
        }

        // get rows and columns of matrix
        struct matrix_slice_left {
            matrix_slice_left(const std::size_t idx, const matrix<T>& _m) : idx_(idx), m(_m) {}
            T& operator[](const std::size_t i) const { return m(idx_, i); }
            //T& operator[](const std::size_t i) { return m(idx_, i); }

            const T* begin() { return &m(idx_,0); }
            const T* end() { return &m(idx_,m.dim2()-1); }

            private:
            const std::size_t idx_;
            const matrix<T>& m;
        };

        matrix_slice_left slice_left(const std::size_t dim1) const { return matrix_slice_left(dim1, *this); }

        struct matrix_slice_right {
            matrix_slice_right(const std::size_t idx, const matrix<T>& _m) : idx_(idx), m(_m) {}
            T operator[](const std::size_t i) const { return m(i, idx_); }
            //T& operator[](const std::size_t i) { return m(i, idx_); }

            struct strided_iterator {
                strided_iterator(T* _x, const std::size_t _stride) : x(_x), stride(_stride) { assert(stride > 0); }
                T operator*() const { return *x; }
                //T& operator*() { return *x; }
                strided_iterator operator+(const std::size_t n) { return strided_iterator(x + stride*n, stride); }
                private:
                T* x;
                const std::size_t stride;
            };

            strided_iterator begin() { return strided_iterator(&m(idx_, 0), m.dim2()); }
            strided_iterator end() { return  strided_iterator(&m(idx_, m.dim2()), m.dim2()); }

            private:
            const std::size_t idx_;
            const matrix<T>& m;
        };

        matrix_slice_right slice_right(const std::size_t dim2) const { return matrix_slice_right({dim2, *this}); }

        // should these functions be members?
        // minima along second dimension
        // should be slower than min2, but possibly it is not, because reduce_min is still a fast operation and not the bottleneck. Measure!
        vector<T> min1() const
        {
            // this function only works if T is REAL!
            static_assert(std::is_same<T,REAL>::value, "");
            vector<T> min(dim1());
            if(std::is_same<T,float>::value || std::is_same<T,double>::value) {
                for(std::size_t x1=0; x1<dim1(); ++x1) {
                    min[x1] = col_min(x1);
                }
            } else {
                assert(false);
            }

            return min;
        }

        //minimum along second dimension, such that to each row of matrix v is added.
        vector<T> min1(const vector<T>& v) const
        {
            assert(v.size() == dim2());
            static_assert(std::is_same<T,REAL>::value, "");
            vector<T> min(dim1());
            for(std::size_t x1=0; x1<dim1(); ++x1) {
                min[x1] = col_min(x1, v);
            }
            return min; 
        }

        // minima along first dimension
        vector<T> min2() const
        {
            static_assert(std::is_same<T,REAL>::value, "");
            vector<T> min(dim2());
            std::fill(min.begin(), min.end(), std::numeric_limits<REAL>::infinity());
            for(size_t i=0; i<dim1(); ++i)
                for(size_t j=0; j<dim2(); ++j)
                    min[j] = std::min(min[j], (*this)(i,j));

            return min;
        }

        // possibly make free function!
        vector<T> min2(const vector<T>& v) const
        {
            static_assert(std::is_same<T,REAL>::value, "");
            assert(v.size() == dim1());
            vector<T> min(dim2());
            for(size_t j=0; j<min.size(); ++j)
                min[j] = std::numeric_limits<REAL>::infinity();

            for(size_t i=0; i<dim1(); ++i)
                for(size_t j=0; j<dim2(); ++j)
                    min[j] = std::min(min[j], (*this)(i,j) + v[i]);

            return min;
        }

        T min() const
        {
            return vec_.min();
        }

        T col_min(const std::size_t x1) const
        {
            assert(x1<dim1());
            REAL min = std::numeric_limits<REAL>::infinity();
            for(size_t j=0; j<dim2(); ++j)
                min = std::min((*this)(x1,j), min);
            return min;
        }

        T col_min(const std::size_t x1, const vector<T>& v) const
        {
            assert(x1<dim1());
            REAL min = std::numeric_limits<REAL>::infinity();
            for(size_t j=0; j<dim2(); ++j)
                min = std::min((*this)(x1,j) + v[j], min);
            return min;
        }

    protected:
        vector<T> vec_;
        std::size_t dim2_;
};

template<typename T=REAL>
class tensor3 : public vector<T> { // do zrobienia: remove
    public:
        tensor3(const std::size_t d1, const std::size_t d2, const std::size_t d3) : vector<T>(d1*d2*d3), dim2_(d2), dim3_(d3) {}
        tensor3(const std::size_t d1, const std::size_t d2, const std::size_t d3, const T val) : tensor3<T>(d1,d2,d3) {
            std::fill(this->begin(), this->end(), val);
        }
        tensor3(const tensor3<T>& o) :
            vector<T>(o),
            dim2_(o.dim2_),
            dim3_(o.dim3_)
    {}
        tensor3(tensor3&& o) :
            vector<T>(std::move(o)),
            dim2_(o.dim2_),
            dim3_(o.dim3_)
    {}
        void operator=(const tensor3<T>& o) {
            assert(this->size() == o.size() && o.dim2_ == dim2_ && o.dim3_ == dim3_);
            for(std::size_t i=0; i<o.size(); ++i) { 
                (*this)[i] = o[i]; 
            }
        }
        tensor3& operator+=(const tensor3<T>& o)
        {
            static_assert(std::is_same<T,double>::value,"");
            assert(dim1() == o.dim1());
            assert(dim2() == o.dim2());
            assert(dim3() == o.dim3());
            this->operator+=(o);
            return *this;
        } 
        T& operator()(const std::size_t x1, const std::size_t x2, const std::size_t x3) { 
            assert(x1<dim1() && x2<dim2() && x3<dim3());
            return (*this)[x1*dim2_*dim3_ + x2*dim3_ + x3]; 
        }
        T operator()(const std::size_t x1, const std::size_t x2, const std::size_t x3) const { return (*this)[x1*dim2_*dim3_ + x2*dim3_ + x3]; }
        const std::size_t dim1() const { return this->size()/(dim2_*dim3_); }
        const std::size_t dim2() const { return dim2_; }
        const std::size_t dim3() const { return dim3_; }
    protected:
        const std::size_t dim2_, dim3_;
};

// given fixed first index, the matrices can be of differing sizes
template<typename T>
class tensor3_variable : public vector<T> {
    public:
        template<typename ITERATOR>
            tensor3_variable(ITERATOR size_begin, ITERATOR size_end)
            : dimensions_(std::distance(size_begin, size_end)),
            vector<T>(size(size_begin, size_end), 0)
    {
        std::size_t offset = 0;
        for(std::size_t i=0; i<std::distance(size_begin, size_end); ++i) {
            dimensions_[i][0] = size_begin[i][0];
            dimensions_[i][1] = size_begin[i][1];
            dimensions_[i].offset = offset;
            offset += dimensions_[i][0] * dimensions_[i][1];
        }
    }
        T& operator()(const std::size_t x1, const std::size_t x2, const std::size_t x3) { 
            assert(x1<dim1() && x2<dim2(x1) && x3<dim3(x1));
            const auto offset = dimensions_[x1].offset;
            return (*this)[offset + x2*dim3(x1) + x3]; 
        }
        T operator()(const std::size_t x1, const std::size_t x2, const std::size_t x3) const 
        {
            assert(x1<dim1() && x2<dim2(x1) && x3<dim3(x1));
            const auto offset = dimensions_[x1].offset;
            return (*this)[offset + x2*dim3(x1) + x3]; 
        }      

        const std::size_t dim1() const { return dimensions_.size(); }
        const std::size_t dim2(const std::size_t x1) const
        {
            assert(x1<dim1());
            return dimensions_[x1][0];
        }
        const std::size_t dim3(const std::size_t x1) const
        {
            assert(x1<dim1());
            return dimensions_[x1][1];
        }

        vector<T>& data() { return *static_cast<vector<T>*>(this); }
        const vector<T>& data() const { return *static_cast<vector<T>*>(this); }

    private:
        template<typename ITERATOR>
            static std::size_t size(ITERATOR size_begin, ITERATOR size_end)
            {
                std::size_t s=0;
                for(auto size_it=size_begin; size_it!=size_end; ++size_it) {
                    s += (*size_it)[0] * (*size_it)[1];
                }
                return s;
            }
        struct dim_entry: public std::array<std::size_t,2> { // dim2, dim3
            std::size_t offset; // offset into data array
        };

        vector<dim_entry> dimensions_; 

};

template<std::size_t FIXED_DIM, typename T=REAL>
class matrix_view_of_tensor : public vector_expression<T,matrix_view_of_tensor<FIXED_DIM,T>> {
    public:
        matrix_view_of_tensor(tensor3<T>& t, const std::size_t fixed_index) : fixed_index_(fixed_index), t_(t) {}
        ~matrix_view_of_tensor() {
            static_assert(FIXED_DIM < 3,"");
        }
        const std::size_t size() const { 
            if(FIXED_DIM==0) {
                return t_.dim2()*t_.dim3();
            } else if(FIXED_DIM == 1) {
                return t_.dim1()*t_.dim2();
            } else {
                return t_.dim2()*t_.dim3();
            }
        }
        const std::size_t dim1() const { 
            if(FIXED_DIM==0) {
                return t_.dim2();
            } else {
                return t_.dim1();
            }
        }
        const std::size_t dim2() const { 
            if(FIXED_DIM==2) {
                return t_.dim2();
            } else {
                return t_.dim3();
            }
        }
        T& operator()(const std::size_t x1, const std::size_t x2) { 
            if(FIXED_DIM==0) {
                return t_(fixed_index_,x1,x2);
            } else if(FIXED_DIM == 1) {
                return t_(x1,fixed_index_,x2);
            } else {
                return t_(x1,x2,fixed_index_);
            }
        }
        T operator()(const std::size_t x1, const std::size_t x2) const { 
            if(FIXED_DIM==0) {
                return t_(fixed_index_,x1,x2);
            } else if(FIXED_DIM == 1) {
                return t_(x1,fixed_index_,x2);
            } else {
                return t_(x1,x2,fixed_index_);
            }
        }

    private:
        const std::size_t fixed_index_;
        tensor3<T>& t_;
};

// primitive expression templates for all the above linear algebraic classes
template<typename T, typename E>
struct scaled_vector : public vector_expression<T,scaled_vector<T,E>> {
    scaled_vector(const T& omega, const E& a) : omega_(omega), a_(a) {}
    const T operator[](const std::size_t i) const {
        return omega_*a_[i];
    }
    const T operator()(const std::size_t i, const std::size_t j) const {
        return omega_*a_(i,j);
    }
    const T operator()(const std::size_t i, const std::size_t j, const std::size_t k) const {
        return omega_*a_(i,j,k);
    }
    std::size_t size() const { return a_.size(); }
    std::size_t dim1() const { return a_.dim1(); }
    std::size_t dim2() const { return a_.dim2(); }
    std::size_t dim3() const { return a_.dim3(); }
    private:
    const T omega_;
    const E& a_;
};


template<typename T, typename E>
scaled_vector<T,vector_expression<T,E>> 
operator*(const T omega, const vector_expression<T,E> & v) {
    return scaled_vector<T,vector_expression<T,E>>(omega, v);
}

template<typename T, typename E>
struct scaled_matrix : public matrix_expression<T,scaled_matrix<T,E>> {
    scaled_matrix(const T& omega, const E& a) : omega_(omega), a_(a) {}
    const T operator[](const std::size_t i) const {
        return omega_*a_[i];
    }
    const T operator()(const std::size_t i, const std::size_t j) const {
        return omega_*a_(i,j);
    }
    const T operator()(const std::size_t i, const std::size_t j, const std::size_t k) const {
        return omega_*a_(i,j,k);
    }
    std::size_t size() const { return a_.size(); }
    std::size_t dim1() const { return a_.dim1(); }
    std::size_t dim2() const { return a_.dim2(); }
    std::size_t dim3() const { return a_.dim3(); }
    private:
    const T omega_;
    const E& a_;
};

template<typename T, typename E>
scaled_matrix<T,matrix_expression<T,E>> 
operator*(const T omega, const matrix_expression<T,E> & v) {
    return scaled_matrix<T,matrix_expression<T,E>>(omega, v);
}


template<typename T, typename E>
struct minus_vector : public vector_expression<T,minus_vector<T,E>> {
    minus_vector(const E& a) : a_(a) {}
    const T operator[](const std::size_t i) const {
        return -a_[i];
    }
    const T operator()(const std::size_t i, const std::size_t j) const {
        return -a_(i,j);
    }
    const T operator()(const std::size_t i, const std::size_t j, const std::size_t k) const {
        return -a_(i,j,k);
    }
    std::size_t size() const { return a_.size(); }
    std::size_t dim1() const { return a_.dim1(); }
    std::size_t dim2() const { return a_.dim2(); }
    std::size_t dim3() const { return a_.dim3(); }
    private:
    const E& a_;
};

template<typename T, typename E>
struct minus_matrix : public matrix_expression<T,minus_matrix<T,E>> {
    minus_matrix(const E& a) : a_(a) {}
    const T operator[](const std::size_t i) const {
        return -a_[i];
    }
    const T operator()(const std::size_t i, const std::size_t j) const {
        return -a_(i,j);
    }
    std::size_t size() const { return a_.size(); }
    std::size_t dim1() const { return a_.dim1(); }
    std::size_t dim2() const { return a_.dim2(); }
    private:
    const E& a_;
};

template<typename T, typename E>
minus_vector<T,vector_expression<T,E>> 
operator-(vector_expression<T,E> const& v) {
    return minus_vector<T,vector_expression<T,E>>(v);
}

template<typename T, typename E>
minus_matrix<T,matrix_expression<T,E>> 
operator-(matrix_expression<T,E> const& v) {
    return minus_matrix<T,matrix_expression<T,E>>(v);
}


} // end namespace LPMP

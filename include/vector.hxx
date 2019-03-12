#ifndef LPMP_VECTOR_HXX
#define LPMP_VECTOR_HXX

#include "memory_allocator.hxx"
//#include "serialization.hxx"
#include "config.hxx"
#include "help_functions.hxx"
//#include "cereal/archives/binary.hpp"

namespace LPMP {

// fixed size vector allocated from block allocator
// possibly holding size explicitly is not needed: It is held by allocator as well

// primitive expression templates for vector
template<typename T, typename E>
class vector_expression {
public:
   T operator[](const INDEX i) const { return static_cast<E const&>(*this)[i]; }
   T operator()(const INDEX i1, const INDEX i2) const { return static_cast<E const&>(*this)(i1,i2); }
   T operator()(const INDEX i1, const INDEX i2, const INDEX i3) const { return static_cast<E const&>(*this)(i1,i2,i3); }
   INDEX size() const { return static_cast<E const&>(*this).size(); }
   INDEX dim1() const { return static_cast<E const&>(*this).dim1(); }
   INDEX dim2() const { return static_cast<E const&>(*this).dim2(); }
   INDEX dim3() const { return static_cast<E const&>(*this).dim3(); }

   E& operator()() { return static_cast<E&>(*this); }
   const E& operator()() const { return static_cast<const E&>(*this); }

   friend std::ostream& operator<<(std::ostream& os, const vector_expression<T,E>& v) {
     for(INDEX i=0; i<v.size(); ++i) {
       os << v[i] << " ";
     }
     os << "\n";
     return os;
   }
};

template<typename T, typename E>
class matrix_expression {
public:
   T operator[](const INDEX i) const { return static_cast<E const&>(*this)[i]; }
   INDEX size() const { return static_cast<E const&>(*this).size(); }

   T operator()(const INDEX i1, const INDEX i2) const { return static_cast<E const&>(*this)(i1,i2); }
   INDEX dim1() const { return static_cast<E const&>(*this).dim1(); }
   INDEX dim2() const { return static_cast<E const&>(*this).dim2(); }

   E& operator()() { return static_cast<E&>(*this); }
   const E& operator()() const { return static_cast<const E&>(*this); }

   friend std::ostream& operator<<(std::ostream& os, const matrix_expression<T,E>& m) {
     for(INDEX i=0; i<m.dim1(); ++i) {
       for(INDEX j=0; j<m.dim2(); ++j) {
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
   T operator[](const INDEX i) const { return static_cast<E const&>(*this)[i]; }
   INDEX size() const { return static_cast<E const&>(*this).size(); }

   T operator()(const INDEX i1, const INDEX i2) const { return static_cast<E const&>(*this)(i1,i2); }
   T operator()(const INDEX i1, const INDEX i2, const INDEX i3) const { return static_cast<E const&>(*this)(i1,i2,i3); }
   INDEX dim1() const { return static_cast<E const&>(*this).dim1(); }
   INDEX dim2() const { return static_cast<E const&>(*this).dim2(); }
   INDEX dim3() const { return static_cast<E const&>(*this).dim3(); }

   E& operator()() { return static_cast<E&>(*this); }
   const E& operator()() const { return static_cast<const E&>(*this); }
};


// possibly also support different allocators: a pure stack allocator without block might be a good choice as well for short-lived memory
template<typename T=REAL>
class vector : public vector_expression<T,vector<T>> {
public:
  template<typename ITERATOR>
  vector(ITERATOR begin, ITERATOR end)
  {
    const std::size_t size = std::distance(begin,end);
    assert(size > 0);
    const std::size_t padding = std::is_same<REAL,T>::value ? (REAL_ALIGNMENT-(size%REAL_ALIGNMENT))%REAL_ALIGNMENT : 0;
    const std::size_t alignment = std::is_same<REAL,T>::value ? 32 : 0;
    //begin_ = (T*) global_real_block_allocator_array[stack_allocator_index].allocate(size+padding,32);
    begin_ = (T*) global_real_block_arena_array[stack_allocator_index].allocate((size+padding)*sizeof(T),32);
    assert(begin_ != nullptr);
    end_ = begin_ + size;
    for(auto it=this->begin(); begin!=end; ++begin, ++it) {
      (*it) = *begin;
    }
    if(padding != 0) {
      std::fill(end_, end_ + padding, std::numeric_limits<REAL>::infinity());
    }
  }

  vector(const INDEX size) 
  {
    const INDEX padding = std::is_same<REAL,T>::value ? (REAL_ALIGNMENT-(size%REAL_ALIGNMENT))%REAL_ALIGNMENT : 0;
    if(std::is_same<REAL,T>::value) {
       assert((padding + size)%REAL_ALIGNMENT == 0);
    }
    //begin_ = global_real_block_allocator_array[stack_allocator_index].allocate(size+padding,32);
    //begin_ = (T*) global_real_block_allocator_array[stack_allocator_index].allocate(size+padding,32);
    begin_ = (T*) global_real_block_arena_array[stack_allocator_index].allocate((size+padding)*sizeof(T),32);
    assert(size > 0);
    assert(begin_ != nullptr);
    end_ = begin_ + size;
    // infinities in padding
    pad_infinity<T>();
  }
  vector(const INDEX size, const T value) 
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
        //global_real_block_allocator_array[stack_allocator_index].deallocate((void*)begin_,1);
        global_real_block_arena_array[stack_allocator_index].deallocate((void*)begin_);
     }
     static_assert(sizeof(T) % sizeof(int) == 0,"");
  }
   vector(const vector& o)  
      : vector(o.size())
   {
      //const INDEX size = o.size();
      //const INDEX padding = std::is_same<REAL,T>::value ? (REAL_ALIGNMENT-(size%REAL_ALIGNMENT))%REAL_ALIGNMENT : 0;
     //begin_ = global_real_block_allocator_array[stack_allocator_index].allocate(o.size()+padding,32);
     //begin_ = (T*) global_real_block_allocator_array[stack_allocator_index].allocate(size+padding,32);
     //begin_ = (T*) global_real_block_arena_array[stack_allocator_index].allocate((size+padding)*sizeof(T),32);
     //end_ = begin_ + o.size();
     assert(begin_ != nullptr);
     auto it = begin_;
     for(auto o_it = o.begin(); o_it!=o.end(); ++it, ++o_it) { *it = *o_it; }
     //pad_infinity<T>();
   }
   template<typename Q=T>
   typename std::enable_if<std::is_same<Q,REAL>::value>::type pad_infinity()
   {
      const INDEX padding = (REAL_ALIGNMENT-(size()%REAL_ALIGNMENT))%REAL_ALIGNMENT;
      std::fill(end_, end_+padding, std::numeric_limits<REAL>::infinity());
   }
   template<typename Q=T>
   typename std::enable_if<!std::is_same<Q,REAL>::value>::type pad_infinity()
   {}

   vector(vector&& o)
   {
      begin_ = o.begin_;
      end_ = o.end_;
      o.begin_ = nullptr;
      o.end_ = nullptr;
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
      for(INDEX i=0; i<size(); i+=REAL_ALIGNMENT) {
         REAL_VECTOR tmp = simdpp::load( begin_+i );
         REAL_VECTOR v = simdpp::load(o.begin() + i);
         simdpp::store(begin_ + i, tmp + v);
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
      for(INDEX i=0; i<o.size(); ++i) { 
         if((*this)[i] != o[i]) { 
            return false;
         }
      }
      return true;
   }
   template<typename E>
   void operator=(const vector_expression<T,E>& o) {
      assert(size() == o.size());
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] = o[i]; 
      }
   }
   template<typename E>
   void operator-=(const vector_expression<T,E>& o) {
      assert(size() == o.size());
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] -= o[i]; } 
   }
   template<typename E>
   void operator+=(const vector_expression<T,E>& o) {
      assert(size() == o.size());
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] += o[i]; } 
   }

   void operator+=(const REAL x) {
      for(INDEX i=0; i<this->size(); ++i) { 
         (*this)[i] += x;
      } 
   }

   void operator-=(const REAL x) {
      for(INDEX i=0; i<this->size(); ++i) { 
         (*this)[i] -= x;
      } 
   }


   // force construction from expression template
   template<typename E>
   vector(vector_expression<T,E>& v) : vector(v.size())
   {
      for(INDEX i=0; v.size(); ++i) {
         (*this)[i] = v[i];
      }
   }

   INDEX size() const { return end_ - begin_; }

   const T& operator[](const INDEX i) const {
      assert(i<size());
      //assert(!std::isnan(begin_[i]));
      return begin_[i];
   }
   T& operator[](const INDEX i) {
      assert(i<size());
      //assert(!std::isnan(begin_[i]));
      return begin_[i];
   }
   using iterator = T*;
   T* begin() const { return begin_; }
   T* end() const { return end_; }

   T back() const { return *(end_-1); }
   T& back() { return *(end_-1); }

   template<typename ARCHIVE>
   void serialize(ARCHIVE& ar)
   {
      assert(false);
      //ar( cereal::binary_data( begin_, sizeof(T)*size()) );
      //ar( binary_data<REAL>( begin_, sizeof(T)*size()) );
   }

   void prefetch() const { simdpp::prefetch_read(begin_); }

   // minimum operation with simd instructions (when T is float, double or integer)
   T min() const
   {
     //check for correct alignment
     //std::cout << std::size_t(begin_) << " aligned?" << std::endl;
     //if((std::size_t(begin_) % 32) != 0) {
     //  std::cout << std::size_t(begin_) << "not aligned" << std::endl;
     //  assert(false);
     //}

     assert((std::size_t(begin_) % 32) == 0);

     if(std::is_same<T,double>::value) {

       REAL_VECTOR min_val = simdpp::load( begin_ );
       for(auto it=begin_+REAL_ALIGNMENT; it<end_; it+=REAL_ALIGNMENT) {
         REAL_VECTOR tmp = simdpp::load( it );
         min_val = simdpp::min(min_val, tmp); 
       }
       return simdpp::reduce_min(min_val);

     } else {
       assert(false);
       return *std::min_element(begin(), end());
     }
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

   // TODO: remove
   void min(const T val)
   {
     assert(false); // should not be used
     static_assert(std::is_same<T,float>::value || std::is_same<T,double>::value,"");
     if(std::is_same<T,float>::value) {
       simdpp::float32<8> val_vec = simdpp::make_float(val);
       for(auto it=begin_+8; it<end_; it+=8) {
         simdpp::float32<8> tmp = simdpp::load( it );
         tmp = simdpp::min(val_vec, tmp); 
         simdpp::store(it, tmp);
       }
     } else if(std::is_same<T,double>::value) {
       simdpp::float64<4> val_vec = simdpp::make_float(val);
       for(auto it=begin_+4; it<end_; it+=4) {
         simdpp::float64<4> tmp = simdpp::load( it );
         tmp = simdpp::min(val_vec, tmp); 
         simdpp::store(it, tmp);
       }
     } else {
       assert(false);
     } 
   }

   std::array<T,2> two_min() const
   {
       assert(size() >= 2);

       if constexpr(std::is_same<T,float>::value) {

           assert(false); // does not work yet
           simdpp::float32<8> min_val = simdpp::make_float(std::numeric_limits<REAL>::infinity());
           simdpp::float32<8> second_min_val = simdpp::make_float(std::numeric_limits<REAL>::infinity());
           for(auto it=begin_; it<end_; it+=8) {
               simdpp::float32<8> tmp = simdpp::load( it );
               simdpp::float32<8> tmp_min = simdpp::min(tmp, min_val);
               simdpp::float32<8> tmp_max = simdpp::max(tmp, min_val);
               min_val = tmp_min;
               second_min_val = simdpp::min(tmp_max, second_min_val); 
           }

           // second minimum is the minimum of the second minimum in vector min and the minimum in second_min_val
           simdpp::float32<4> x1;
           simdpp::float32<4> y1;
           simdpp::split(min_val, x1, y1);
           simdpp::float32<4> min2 = simdpp::min(x1, y1);
           simdpp::float32<4> max2 = simdpp::max(x1, y1);

           std::array<T,4> min_array; //
           min_array[0] = min2.vec(0); //min2.vec(1), min2.vec(2), min2.vec(3)});
           const auto min3 = two_smallest_elements<T>(min_array.begin(), min_array.end());
           //const REAL min = simdpp::reduce_min(min2);
           const REAL second_min = std::min({simdpp::reduce_min(second_min_val), simdpp::reduce_min(max2), min3[1]});
           assert(second_min >= min3[0]);
           return std::array<T,2>({min3[0], second_min});

   } else if constexpr(std::is_same<T,double>::value) {

           simdpp::float64<4> min_val = simdpp::make_float(std::numeric_limits<REAL>::infinity());
           simdpp::float64<4> second_min_val = simdpp::make_float(std::numeric_limits<REAL>::infinity());
           for(auto it=begin_; it<end_; it+=4) {
               simdpp::float64<4> tmp = simdpp::load( it );
               simdpp::float64<4> tmp_min = simdpp::min(tmp, min_val);
               simdpp::float64<4> tmp_max = simdpp::max(tmp, min_val);
               min_val = tmp_min;
               second_min_val = simdpp::min(tmp_max, second_min_val); 
           }

           // second minimum is the minimu of the second minimum in vector min and the minimum in second_min_val

           T smallest = std::numeric_limits<T>::infinity();
           T second_smallest = std::numeric_limits<T>::infinity();
           //for(std::size_t i=0; i<4; ++i) {
               {
               const T min = std::min(smallest, simdpp::extract<0>(min_val));
               const T max = std::max(smallest, simdpp::extract<0>(min_val));
               smallest = min;
               second_smallest = std::min(max, second_smallest);
               }
               {
               const T min = std::min(smallest, simdpp::extract<1>(min_val));
               const T max = std::max(smallest, simdpp::extract<1>(min_val));
               smallest = min;
               second_smallest = std::min(max, second_smallest);
               }
               { 
               const T min = std::min(smallest, simdpp::extract<2>(min_val));
               const T max = std::max(smallest, simdpp::extract<2>(min_val));
               smallest = min;
               second_smallest = std::min(max, second_smallest);
               }
               { 
               const T min = std::min(smallest, simdpp::extract<3>(min_val));
               const T max = std::max(smallest, simdpp::extract<3>(min_val));
               smallest = min;
               second_smallest = std::min(max, second_smallest);
               }
           //}
           second_smallest = std::min(second_smallest, simdpp::reduce_min(second_min_val));

           return {smallest, second_smallest};

   } else {
       assert(false);
   }
           /*
              simdpp::float32<2> x2;
              simdpp::float32<2> y2;
              simdpp::split(min2, x2, y2);
     simdpp::float32<2> min3 = simdpp::min(x2, y2);
     simdpp::float32<2> max3 = simdpp::max(x2, y2);


     const REAL min = simdpp::reduce_min(min3);
     const REAL second_min = std::min({simdpp::reduce_min(second_min_val), simdpp::reduce_min(max2), simdpp::reduce_min(max3), simdpp::reduce_max(min3)});
     assert(min == two_smallest_elements<T>(begin_, end_)[0]);
     assert(second_min == two_smallest_elements<T>(begin_, end_)[1]);
     return std::array<T,2>({min, second_min});
     */
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

    small_vector(const INDEX size)
    {
        allocate(size);
    }

    small_vector(const INDEX size, const T val)
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
    alignas(REAL_ALIGNMENT*sizeof(REAL)) std::array<T,N> local_;
};


template<typename T, INDEX N>
using array_impl = std::array<T,N>;



template<typename T, INDEX N>
class array : public vector_expression<T,array<T,N>> {
public:
   array() {}

   template<typename ITERATOR>
   array(ITERATOR begin, ITERATOR end)
   {
      const INDEX size = std::distance(begin,end);
      assert(size == N);
      for(auto it=array_.begin(); it!=array_.end(); ++it) {
         (*it) = *begin;
      }
   }

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
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] -= o[i]; } 
   }
   template<typename E>
   void operator+=(const vector_expression<T,E>& o) {
      assert(size() == o.size());
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] += o[i]; } 
   }


   // force construction from expression template
   template<typename E>
   array(vector_expression<T,E>& v) 
   {
      for(INDEX i=0; v.size(); ++i) {
         (*this)[i] = v[i];
      }
   }

   constexpr static std::size_t size() { return N; }

   T operator[](const INDEX i) const {
      assert(i<size());
      //assert(!std::isnan(array_[i]));
      return array_[i];
   }
   T& operator[](const INDEX i) {
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

      if(array_.size() > REAL_ALIGNMENT) {
         REAL_VECTOR cur_min = simdpp::load(&array_[0]);
         INDEX last_aligned = array_.size() - (array_.size()%REAL_ALIGNMENT);
         for(auto i=REAL_ALIGNMENT; i<last_aligned; i+=REAL_ALIGNMENT) {
            const REAL_VECTOR tmp = simdpp::load( &array_[i] );
            cur_min = simdpp::min(cur_min, tmp); 
         }

         REAL aligned_min = simdpp::reduce_min(cur_min);

         for(INDEX i=last_aligned; i<array_.size(); ++i) {
            aligned_min = std::min(aligned_min, array_[i]);
         }
         assert(aligned_min == *std::min_element(array_.begin(), array_.end()) );
         return aligned_min; 
      } else { // small vector
         return *std::min_element(array_.begin(), array_.end());
      } 
   }
private:
   alignas(REAL_ALIGNMENT*sizeof(REAL)) std::array<T,N> array_; // should only be done for REALs, not generally
};


// matrix is based on vector.
// However, the entries are padded so that each coordinate (i,0) is aligned
template<typename T=REAL>
class matrix : public matrix_expression<T,matrix<T>> {
public:
   T& operator[](const INDEX i) 
   { 
      assert(i<size());
      const std::size_t x1 = i/dim2();
      const std::size_t x2 = i%dim2();
      return (*this)(x1,x2);
   }
   const T operator[](const INDEX i) const 
   { 
      assert(i<size());
      const std::size_t x1 = i/dim2();
      const std::size_t x2 = i%dim2();
      return (*this)(x1,x2);
   }
   INDEX size() const { return dim1()*dim2(); }
   template<typename ARCHIVE>
   void serialize(ARCHIVE& ar)
   {
      ar( vec_ );
   }

   class const_iterator {
   public:
      const_iterator(const INDEX _dim2, const INDEX _padded_dim2, const vector<T>& _vec, const INDEX _i)
         :
            i(_i),
            dim2(_dim2),
            padded_dim2(_padded_dim2),
            vec(_vec)
      {}

      void operator++() {
         ++i; 
      }
      T operator*() const {
         const INDEX x1 = i/dim2;
         const INDEX x2 = i%dim2;
         return vec[x1*padded_dim2 + x2]; 
      }
      bool operator==(const const_iterator& o) const {
         assert(dim2 == o.dim2 && padded_dim2 == o.padded_dim2 && &vec == &o.vec);
         return i == o.i;
      }
      bool operator!=(const const_iterator& o) const {
         return !( (*this) == o );
      }
   private:
      INDEX i;
      const INDEX dim2;
      const INDEX padded_dim2; 
      const vector<T>& vec;
   };

   class iterator {
   public:
      iterator(const INDEX _dim2, const INDEX _padded_dim2, vector<T>& _vec, const INDEX _i)
         :
            i(_i),
            dim2(_dim2),
            padded_dim2(_padded_dim2),
            vec(_vec)
      {}

      void operator++() {
         ++i; 
      }
      T& operator*() {
         const INDEX x1 = i/dim2;
         const INDEX x2 = i%dim2;
         return vec[x1*padded_dim2 + x2];
      }
      T operator*() const {
         const INDEX x1 = i/dim2;
         const INDEX x2 = i%dim2;
         return vec[x1*padded_dim2 + x2]; 
      }
      bool operator==(const iterator& o) const {
         assert(dim2 == o.dim2 && padded_dim2 == o.padded_dim2 && &vec == &o.vec);
         return i == o.i;
      }
      bool operator!=(const iterator& o) const {
         return !( (*this) == o );
      }
   private:
      INDEX i;
      const INDEX dim2;
      const INDEX padded_dim2; 
      vector<T>& vec;
   };

   // when infinity padding is on, we must jump over those places
   const_iterator begin() const { return const_iterator(dim2(), padded_dim2(), vec_, 0); }
   const_iterator end() const { return const_iterator(dim2(), padded_dim2(), vec_, size()); }

   iterator begin() { return iterator(dim2(), padded_dim2(), vec_, 0); }
   iterator end() { return iterator(dim2(), padded_dim2(), vec_, size()); }

   static INDEX underlying_vec_size(const INDEX d1, const INDEX d2)
   {
     return d1*(d2 + padding(d2));
   }
   static INDEX padding(const INDEX i) {
      if constexpr(std::is_same<T,float>::value || std::is_same<T,double>::value,"")
         return (REAL_ALIGNMENT-(i%REAL_ALIGNMENT))%REAL_ALIGNMENT;
      else
         return 0;
   }

   INDEX padded_dim2() const { return padded_dim2_; }

   void fill_padding()
   {
      if constexpr(std::is_same<T,float>::value || std::is_same<T,double>::value,"") {
         for(INDEX x1=0; x1<dim1(); ++x1) {
            for(INDEX x2=dim2(); x2<padded_dim2(); ++x2) {
               vec_[x1*padded_dim2() + x2] = std::numeric_limits<T>::infinity();
            }
         }
      }
   }
   matrix(const INDEX d1, const INDEX d2) : vec_(underlying_vec_size(d1,d2)), dim2_(d2), padded_dim2_(d2 + padding(d2)) {
      assert(d1 > 0 && d2 > 0);
      fill_padding();
   }
   matrix(const INDEX d1, const INDEX d2, const T val) : vec_(underlying_vec_size(d1,d2)), dim2_(d2), padded_dim2_(d2 + padding(d2)) {
      assert(d1 > 0 && d2 > 0);
      std::fill(vec_.begin(), vec_.end(), val);
      fill_padding();
   }
   matrix() : vec_(), dim2_(0), padded_dim2_(0) {}
   matrix(const matrix& o) 
      : vec_(o.vec_),
      dim2_(o.dim2_),
      padded_dim2_(o.padded_dim2_)
   {}
   matrix(matrix&& o) 
      : vec_(std::move(o.vec_)),
      dim2_(o.dim2_),
      padded_dim2_(o.padded_dim2_)
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
       std::swap(padded_dim2_, copy.padded_dim2_);
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
      for(INDEX i=0; i<o.dim1(); ++i) { 
         for(INDEX j=0; j<o.dim2(); ++j) { 
            (*this)(i,j) -= o(i,j); 
         } 
      }
   }


   T& operator()(const INDEX x1, const INDEX x2) { assert(x1<dim1() && x2<dim2()); return vec_[x1*padded_dim2() + x2]; }
   const T& operator()(const INDEX x1, const INDEX x2) const { assert(x1<dim1() && x2<dim2()); return vec_[x1*padded_dim2() + x2]; }
   T& operator()(const INDEX x1, const INDEX x2, const INDEX x3) { assert(x3 == 0); return (*this)(x1,x2); } // sometimes we treat a matrix as a tensor with trivial last dimension
   const T& operator()(const INDEX x1, const INDEX x2, const INDEX x3) const { assert(x3 == 0); return (*this)(x1,x2); } // sometimes we treat a matrix as a tensor with trivial last dimension

   const INDEX dim1() const { return vec_.size()/padded_dim2(); }
   const INDEX dim2() const { return dim2_; }

   void transpose() {
      assert(dim1() == dim2());
      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x2=0; x2<x1; ++x2) {
            std::swap((*this)(x1,x2), (*this)(x2,x1));
         }
      }
   }

   // get rows and columns of matrix
   struct matrix_slice_left {
      matrix_slice_left(const INDEX idx, const matrix<T>& _m) : idx_(idx), m(_m) {}
      T& operator[](const INDEX i) const { return m(idx_, i); }
      //T& operator[](const INDEX i) { return m(idx_, i); }

      const T* begin() { return &m(idx_,0); }
      const T* end() { return &m(idx_,m.dim2()-1); }

      private:
      const INDEX idx_;
      const matrix<T>& m;
   };

   matrix_slice_left slice_left(const INDEX dim1) const { return matrix_slice_left(dim1, *this); }

   struct matrix_slice_right {
      matrix_slice_right(const INDEX idx, const matrix<T>& _m) : idx_(idx), m(_m) {}
      T operator[](const INDEX i) const { return m(i, idx_); }
      //T& operator[](const INDEX i) { return m(i, idx_); }

      struct strided_iterator {
         strided_iterator(T* _x, const INDEX _stride) : x(_x), stride(_stride) { assert(stride > 0); }
         T operator*() const { return *x; }
         //T& operator*() { return *x; }
         strided_iterator operator+(const INDEX n) { return strided_iterator(x + stride*n, stride); }
         private:
         T* x;
         const INDEX stride;
      };

      strided_iterator begin() { return strided_iterator(&m(idx_, 0), m.padded_dim2()); }
      strided_iterator end() { return  strided_iterator(&m(idx_, m.padded_dim2()), m.padded_dim2()); }

      private:
      const INDEX idx_;
      const matrix<T>& m;
   };

   matrix_slice_right slice_right(const INDEX dim2) const { return matrix_slice_right({dim2, *this}); }

   // should these functions be members?
   // minima along second dimension
   // should be slower than min2, but possibly it is not, because reduce_min is still a fast operation and not the bottleneck. Measure!
   vector<T> min1() const
   {
     // this function only works if T is REAL!
     static_assert(std::is_same<T,REAL>::value, "");
     vector<T> min(dim1());
     if(std::is_same<T,float>::value || std::is_same<T,double>::value) {
       for(INDEX x1=0; x1<dim1(); ++x1) {
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
      for(INDEX x1=0; x1<dim1(); ++x1) {
         min[x1] = col_min(x1, v);
      }
      return min; 
   }

   // minima along first dimension
   vector<T> min2() const
   {
     static_assert(std::is_same<T,REAL>::value, "");
     vector<T> min(dim2());
     // possibly iteration strategy is faster, e.g. doing a non-contiguous access, or explicitly holding a few variables and not storing them back in vector min for a few sizes
     if(std::is_same<T,float>::value || std::is_same<T,double>::value) {
       for(INDEX x2=0; x2<dim2(); x2+=REAL_ALIGNMENT) {
         REAL_VECTOR tmp = simdpp::load( vec_.begin() + x2 );
         simdpp::store(&min[x2], tmp);
       }

       for(INDEX x1=1; x1<dim1(); ++x1) {
         for(INDEX x2=0; x2<dim2(); x2+=REAL_ALIGNMENT) {
           REAL_VECTOR tmp = simdpp::load( vec_.begin() + x1*padded_dim2() + x2 );
           REAL_VECTOR cur_min = simdpp::load( &min[x2] );
           auto updated_min = simdpp::min(cur_min, tmp);
           simdpp::store(&min[x2], updated_min);
         } 
       }
     } else {
       assert(false);
     }

     return min; 
   }

   // possibly make free function!
   vector<T> min2(const vector<T>& v) const
   {
     static_assert(std::is_same<T,REAL>::value, "");
     assert(v.size() == dim1());
     vector<T> min(dim2());

     for(INDEX x2=0; x2<dim2(); x2+=REAL_ALIGNMENT) {
        REAL_VECTOR tmp = simdpp::load( vec_.begin() + x2 );
        REAL_VECTOR _v = simdpp::load_splat(v.begin());
        REAL_VECTOR _sum = tmp + _v;
        simdpp::store(&min[x2], _sum);
     }
     for(INDEX x1=1; x1<dim1(); ++x1) {
        for(INDEX x2=0; x2<dim2(); x2+=REAL_ALIGNMENT) {
           REAL_VECTOR tmp = simdpp::load( vec_.begin() + x1*padded_dim2() + x2 );
           REAL_VECTOR _v = simdpp::load_splat(v.begin() + x1);
           REAL_VECTOR _sum = tmp + _v;
           REAL_VECTOR cur_min = simdpp::load( &min[x2] );
           auto updated_min = simdpp::min(cur_min, _sum);
           simdpp::store(&min[x2], updated_min);
        } 
     }

     return min;
   }

   T min() const
   {
     return vec_.min();
   }

   T col_min(const INDEX x1) const
   {
     assert(x1<dim1());
     REAL_VECTOR cur_min = simdpp::load( vec_.begin() + x1*padded_dim2() );
     for(INDEX x2=REAL_ALIGNMENT; x2<dim2(); x2+=REAL_ALIGNMENT) {
       REAL_VECTOR tmp = simdpp::load( vec_.begin() + x1*padded_dim2() + x2 );
       cur_min = simdpp::min(cur_min, tmp); 
     }
     return simdpp::reduce_min(cur_min); 
   }

   T col_min(const INDEX x1, const vector<T>& v) const
   {
     assert(x1<dim1());
     REAL_VECTOR cur_min = simdpp::load( vec_.begin() + x1*padded_dim2() );
     REAL_VECTOR _v = simdpp::load(v.begin());
     cur_min = cur_min + _v;
     for(INDEX x2=REAL_ALIGNMENT; x2<dim2(); x2+=REAL_ALIGNMENT) {
       REAL_VECTOR tmp = simdpp::load( vec_.begin() + x1*padded_dim2() + x2 );
       REAL_VECTOR _v = simdpp::load(v.begin() + x2);
       tmp = tmp + _v;
       cur_min = simdpp::min(cur_min, tmp); 
     }
     return simdpp::reduce_min(cur_min); 
   }

protected:
   vector<T> vec_;
   INDEX dim2_;
   INDEX padded_dim2_; // possibly do not store but compute when needed?
};

template<typename T=REAL>
class tensor3 : public vector<T> { // do zrobienia: remove
public:
   tensor3(const INDEX d1, const INDEX d2, const INDEX d3) : vector<T>(d1*d2*d3), dim2_(d2), dim3_(d3) {}
   tensor3(const INDEX d1, const INDEX d2, const INDEX d3, const T val) : tensor3<T>(d1,d2,d3) {
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
      for(INDEX i=0; i<o.size(); ++i) { 
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
   T& operator()(const INDEX x1, const INDEX x2, const INDEX x3) { 
      assert(x1<dim1() && x2<dim2() && x3<dim3());
      return (*this)[x1*dim2_*dim3_ + x2*dim3_ + x3]; 
   }
   T operator()(const INDEX x1, const INDEX x2, const INDEX x3) const { return (*this)[x1*dim2_*dim3_ + x2*dim3_ + x3]; }
   const INDEX dim1() const { return this->size()/(dim2_*dim3_); }
   const INDEX dim2() const { return dim2_; }
   const INDEX dim3() const { return dim3_; }
protected:
   const INDEX dim2_, dim3_;
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
    T& operator()(const INDEX x1, const INDEX x2, const INDEX x3) { 
        assert(x1<dim1() && x2<dim2(x1) && x3<dim3(x1));
        const auto offset = dimensions_[x1].offset;
        return (*this)[offset + x2*dim3(x1) + x3]; 
    }
    T operator()(const INDEX x1, const INDEX x2, const INDEX x3) const 
    {
        assert(x1<dim1() && x2<dim2(x1) && x3<dim3(x1));
        const auto offset = dimensions_[x1].offset;
        return (*this)[offset + x2*dim3(x1) + x3]; 
    }      

    const INDEX dim1() const { return dimensions_.size(); }
    const INDEX dim2(const INDEX x1) const
    {
        assert(x1<dim1());
        return dimensions_[x1][0];
    }
    const INDEX dim3(const INDEX x1) const
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

template<INDEX FIXED_DIM, typename T=REAL>
class matrix_view_of_tensor : public vector_expression<T,matrix_view_of_tensor<FIXED_DIM,T>> {
public:
   matrix_view_of_tensor(tensor3<T>& t, const INDEX fixed_index) : fixed_index_(fixed_index), t_(t) {}
   ~matrix_view_of_tensor() {
      static_assert(FIXED_DIM < 3,"");
   }
   const INDEX size() const { 
      if(FIXED_DIM==0) {
         return t_.dim2()*t_.dim3();
      } else if(FIXED_DIM == 1) {
         return t_.dim1()*t_.dim2();
      } else {
         return t_.dim2()*t_.dim3();
      }
   }
   const INDEX dim1() const { 
      if(FIXED_DIM==0) {
         return t_.dim2();
      } else {
         return t_.dim1();
      }
   }
   const INDEX dim2() const { 
      if(FIXED_DIM==2) {
         return t_.dim2();
      } else {
         return t_.dim3();
      }
   }
   T& operator()(const INDEX x1, const INDEX x2) { 
      if(FIXED_DIM==0) {
         return t_(fixed_index_,x1,x2);
      } else if(FIXED_DIM == 1) {
         return t_(x1,fixed_index_,x2);
      } else {
         return t_(x1,x2,fixed_index_);
      }
   }
   T operator()(const INDEX x1, const INDEX x2) const { 
      if(FIXED_DIM==0) {
         return t_(fixed_index_,x1,x2);
      } else if(FIXED_DIM == 1) {
         return t_(x1,fixed_index_,x2);
      } else {
         return t_(x1,x2,fixed_index_);
      }
   }

private:
   const INDEX fixed_index_;
   tensor3<T>& t_;
};

// primitive expression templates for all the above linear algebraic classes
template<typename T, typename E>
struct scaled_vector : public vector_expression<T,scaled_vector<T,E>> {
   scaled_vector(const T& omega, const E& a) : omega_(omega), a_(a) {}
   const T operator[](const INDEX i) const {
      return omega_*a_[i];
   }
   const T operator()(const INDEX i, const INDEX j) const {
      return omega_*a_(i,j);
   }
   const T operator()(const INDEX i, const INDEX j, const INDEX k) const {
      return omega_*a_(i,j,k);
   }
   INDEX size() const { return a_.size(); }
   INDEX dim1() const { return a_.dim1(); }
   INDEX dim2() const { return a_.dim2(); }
   INDEX dim3() const { return a_.dim3(); }
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
   const T operator[](const INDEX i) const {
      return omega_*a_[i];
   }
   const T operator()(const INDEX i, const INDEX j) const {
      return omega_*a_(i,j);
   }
   const T operator()(const INDEX i, const INDEX j, const INDEX k) const {
      return omega_*a_(i,j,k);
   }
   INDEX size() const { return a_.size(); }
   INDEX dim1() const { return a_.dim1(); }
   INDEX dim2() const { return a_.dim2(); }
   INDEX dim3() const { return a_.dim3(); }
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
   const T operator[](const INDEX i) const {
      return -a_[i];
   }
   const T operator()(const INDEX i, const INDEX j) const {
      return -a_(i,j);
   }
   const T operator()(const INDEX i, const INDEX j, const INDEX k) const {
      return -a_(i,j,k);
   }
   INDEX size() const { return a_.size(); }
   INDEX dim1() const { return a_.dim1(); }
   INDEX dim2() const { return a_.dim2(); }
   INDEX dim3() const { return a_.dim3(); }
   private:
   const E& a_;
};

template<typename T, typename E>
struct minus_matrix : public matrix_expression<T,minus_matrix<T,E>> {
   minus_matrix(const E& a) : a_(a) {}
   const T operator[](const INDEX i) const {
      return -a_[i];
   }
   const T operator()(const INDEX i, const INDEX j) const {
      return -a_(i,j);
   }
   INDEX size() const { return a_.size(); }
   INDEX dim1() const { return a_.dim1(); }
   INDEX dim2() const { return a_.dim2(); }
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

#endif // LPMP_VECTOR_HXX

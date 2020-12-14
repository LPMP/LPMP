#pragma once

#include <vector>
#include <array>
#include <cassert>

namespace LPMP {

// general two-dimensional array with variable first and second dimension sizes, i.e. like vector<vector<T>>. Holds all data contiguously and therefore may be more efficient than vector<vector<T>>

// do zrobienia: - alignment?
//               - custom memory allocator?
//               - iterators
template<typename T>
class two_dim_variable_array
{
public:
    template<class I> friend class two_dim_variable_array;

   two_dim_variable_array() 
       : two_dim_variable_array(std::vector<std::size_t>{})
   {}

   ~two_dim_variable_array() {
      static_assert(!std::is_same_v<T,bool>, "value type cannot be bool");
   }

   template<typename I>
       two_dim_variable_array(const two_dim_variable_array<I>& o)
       : offsets_(o.offsets_),
       data_(offsets_.back())
       {}

   two_dim_variable_array(const two_dim_variable_array<T>& o)
       : offsets_(o.offsets_),
       data_(o.data_)
    {}

   // iterator holds size of each dimension of the two dimensional array
   template<typename I>
   two_dim_variable_array(const std::vector<I>& size)
   : offsets_(compute_offsets(size.begin(), size.end())),
    data_(offsets_.back())
   {
      //const std::size_t s = set_dimensions(size.begin(), size.end());
      //data_.resize(s);
   }

   two_dim_variable_array(const std::vector<std::vector<T>>& o)
   {
      std::vector<std::size_t> size_array;
      size_array.reserve(o.size());
      for(const auto& v : o) {
         size_array.push_back(v.size());
      }
      const std::size_t s = set_dimensions(size_array.begin(), size_array.end());
      data_.resize(s);

      for(std::size_t i=0; i<size(); ++i) {
         for(std::size_t j=0; j<(*this)[i].size(); ++j) {
            (*this)(i,j) = o[i][j];
         }
      }
   }

   // iterator holds size of each dimension of the two dimensional array
   template<typename ITERATOR>
   two_dim_variable_array(ITERATOR size_begin, ITERATOR size_end)
   : offsets_(compute_offsets(size_begin, size_end)),
    data_(offsets_.back())
   {
      //const std::size_t s = set_dimensions(size_begin, size_end);
      //data_.resize(s);
   }

   // iterator holds size of each dimension of the two dimensional array
   template<typename ITERATOR>
   two_dim_variable_array(ITERATOR size_begin, ITERATOR size_end, const T& val)
   : offsets_(compute_offsets(size_begin, size_end)),
    data_(offsets_.back(), val)
   {
      //const std::size_t s = set_dimensions(size_begin, size_end);
      //data_.resize(s);
   }

   /*
   friend std::ostream& operator<<(std::ostream& os, const two_dim_variable_array<T>& a) {
     for(std::size_t i=0; i<a.size(); ++i) {
       for(std::size_t j=0; j<a[i].size(); ++j) {
         os << a(i,j) << " ";
       }
       os << "\n";
     }
     return os;
   }
   */

   void clear();
   template<typename ITERATOR>
       void push_back(ITERATOR val_begin, ITERATOR val_end);

   template<typename ITERATOR>
   void resize(ITERATOR begin, ITERATOR end)
   {
      const std::size_t size = set_dimensions(begin,end);
      data_.resize(size); 
   }
   template<typename ITERATOR>
   void resize(ITERATOR begin, ITERATOR end, T val)
   {
      const std::size_t size = set_dimensions(begin,end);
      data_.resize(size, val); 
   }

   struct ConstArrayAccessObject
   {
      ConstArrayAccessObject(const T* begin, const T* end) : begin_(begin), end_(end) { assert(begin <= end); }
      const T& operator[](const std::size_t i) const { assert(i < size()); return begin_[i]; }
      std::size_t size() const {  return (end_ - begin_); }

      const T* begin() { return begin_; }
      const T* end() { return end_; }
      auto rbegin() { return std::make_reverse_iterator(end()); }
      auto rend() { return std::make_reverse_iterator(begin()); }
      const T& back() { return *(end()-1); }

      const T* begin() const { return begin_; }
      const T* end() const { return end_; }
      auto rbegin() const { return std::make_reverse_iterator(end()); }
      auto rend() const { return std::make_reverse_iterator(begin()); }
      const T& back() const { return *(end()-1); }

      private:
         const T* begin_;
         const T* end_;
   };
   struct ArrayAccessObject
   {
      ArrayAccessObject(T* begin, T* end) : begin_(begin), end_(end) { assert(begin <= end); }
      template<typename VEC>
      void operator=(const VEC& o) 
      { 
        assert(o.size() == this->size());
        const auto s = this->size();
        for(std::size_t i=0; i<s; ++i) {
          (*this)[i] = o[i];
        }
      }
      const T& operator[](const std::size_t i) const { assert(i < size()); return begin_[i]; }
      T& operator[](const std::size_t i) { assert(i < size()); return begin_[i]; }
      std::size_t size() const {  return (end_ - begin_); }

      T* begin() { return begin_; }
      T* end() { return end_; }
      auto rbegin() { return std::make_reverse_iterator(end()); }
      auto rend() { return std::make_reverse_iterator(begin()); }
      T& back() { return *(end()-1); }

      T const* begin() const { return begin_; }
      T const* end() const { return end_; }
      auto rbegin() const { return std::make_reverse_iterator(end()); }
      auto rend() const { return std::make_reverse_iterator(begin()); }
      const T& back() const { return *(end()-1); }

      private:
         T* begin_;
         T* end_;
   };
   ArrayAccessObject operator[](const std::size_t i) {
      assert(i<this->size());
      return ArrayAccessObject( &data_[offsets_[i]], &data_[offsets_[i+1]] );
   }
   ConstArrayAccessObject operator[](const std::size_t i) const {
      assert(i<this->size());
      return ConstArrayAccessObject( &data_[offsets_[i]], &data_[offsets_[i+1]] );
   }
   const T& operator()(const std::size_t i, const std::size_t j) const
   {
      assert(i<size() && j< (*this)[i].size());
      return data_[offsets_[i]+j];
   }
   T& operator()(const std::size_t i, const std::size_t j)
   {
      assert(i < size() && j< (*this)[i].size());
      return data_[offsets_[i]+j];
   }

   std::size_t no_elements() const { return data_.size(); }
   std::size_t size() const { assert(offsets_.size() > 0); return offsets_.size()-1; }

   ConstArrayAccessObject back() const 
   {
       const std::size_t i = size()-1;
       return (*this)[i];
   }

   ArrayAccessObject back()
   {
       const std::size_t i = size()-1;
       return (*this)[i];
   }

   struct iterator : public std::iterator< std::random_access_iterator_tag, T* > {
     iterator(T* _data, std::size_t* _offset) : data(_data), offset(_offset) {}
     void operator++() { ++offset; }
     void operator--() { --offset; }
     iterator& operator+=(const std::size_t i) { offset+=i; return *this; }
     iterator& operator-=(const std::size_t i) { offset-=i; return *this; }
     iterator operator+(const std::size_t i) { iterator it(data,offset+i); return it; }
     iterator operator-(const std::size_t i) { iterator it(data,offset-i); return it; }
     auto operator-(const iterator it) const { return offset - it.offset; }
     ArrayAccessObject operator*() { return ArrayAccessObject(data+*offset,data+*(offset+1)); }
     const ArrayAccessObject operator*() const { return ArrayAccessObject(data+*offset,data+*(offset+1)); }
     bool operator==(const iterator it) const { return data == it.data && offset == it.offset; }
     bool operator!=(const iterator it) const { return !(*this == it); }
     T* data;
     std::size_t* offset;
   };

   struct reverse_iterator : public std::iterator_traits< T* > {
     reverse_iterator(T* _data, std::size_t* _offset) : data(_data), offset(_offset) {}
     void operator++() { --offset; }
     void operator--() { ++offset; }
     reverse_iterator& operator+=(const std::size_t i) { offset-=i; return *this; }
     reverse_iterator& operator-=(const std::size_t i) { offset+=i; return *this; }
     iterator operator+(const std::size_t i) { iterator it(data,offset-i); return it; }
     iterator operator-(const std::size_t i) { iterator it(data,offset+i); return it; }
     auto operator-(const reverse_iterator it) const { return it.offset - offset; }
     ArrayAccessObject operator*() { return ArrayAccessObject(data+*(offset-1),data+*offset); }
     const ArrayAccessObject operator*() const { return ArrayAccessObject(data+*(offset-1),data+*offset); }
     bool operator==(const reverse_iterator it) const { return data == it.data && offset == it.offset; }
     bool operator!=(const reverse_iterator it) const { return !(*this == it); }
     T* data;
     std::size_t* offset;
   };

   iterator begin() { return iterator(&data_[0],&offsets_[0]); }
   iterator end() { return iterator(&data_[0],&offsets_.back()); }

   // TODO: correct?
   iterator begin() const { return iterator(&data_[0],&offsets_[0]); }
   iterator end() const { return iterator(&data_[0],&offsets_.back()); }

   reverse_iterator rbegin() { return reverse_iterator(&data_[0],&offsets_.back()); }
   reverse_iterator rend() { return reverse_iterator(&data_[0],&offsets_[0]); }

   struct size_iterator : public std::iterator_traits<const std::size_t*> {
      size_iterator(const std::size_t* _offset) : offset(_offset) {}
      std::size_t operator*() const { return *(offset+1) - *offset; } 
      void operator++() { ++offset; }
      auto operator-(const size_iterator it) const { return offset - it.offset; }
      size_iterator operator-(const std::size_t i) const { size_iterator it(offset-i); return it; }
      bool operator==(const size_iterator it) const { return offset == it.offset; }
      bool operator!=(const size_iterator it) const { return !(*this == it); }
      const std::size_t* offset;
   };

   auto size_begin() const { assert(offsets_.size() > 0); return size_iterator({&offsets_[0]}); }
   auto size_end() const { assert(offsets_.size() > 0); return size_iterator({&offsets_.back()}); }

   void set(const T& val)
   {
	   for(std::size_t i=0; i<size(); ++i) {
		   for(std::size_t j=0; j<(*this)[i].size(); ++j) {
			   (*this)(i,j) = val;
		   }
	   }
   }

   auto& data() { return data_; }
   const auto& data() const { return data_; }

   std::size_t first_index(const T* p) const
   {
       assert(false); // TODO: does not work yet, write test
       assert(p >= &data_[0]);
       assert(p < &data_[0] + data_.size());

       // binary search to locate first index
       std::size_t lower = 0;
       std::size_t upper = size()-1;
       while(lower <= upper) {
           const std::size_t middle = (lower+upper)/2;
           if( &(*this)(middle,0) <= p && p <= &(*this)[middle].back()) {
               return middle; 
           } else if( &(*this)(middle,0) < p ) {
               lower = middle+1;
           } else {
               upper = middle-1; 
           }
       }
       assert(false);
       return std::numeric_limits<std::size_t>::max();
   }

   std::array<std::size_t,2> indices(const T* p) const
   {
       assert(p >= &data_[0]);
       assert(p < &data_[0] + data_.size());

       const std::size_t first_index = first_index(p);

       // binary search for second index
       const std::size_t second_index = [&]() {
           std::size_t lower = 0;
           std::size_t upper = (*this)[first_index].size()-1;
           while(lower <= upper) {
               const std::size_t middle = (lower+upper)/2;
               if(&(*this)(first_index,middle) < p)
                   lower = middle+1; 
               else if(&(*this)(first_index,middle) > p)
                   upper = middle-1;
               else
                   return middle;
           }
           assert(false);
           return std::numeric_limits<std::size_t>::max();
       }();

       assert(p == &(*this)(first_index, second_index));
       return {first_index, second_index};
   }

   std::size_t index(const std::size_t first_index, const std::size_t second_index) const
   {
       assert(first_index < size());
       assert(second_index < (*this)[first_index].size());
       return std::distance(&data_[0], &(*this)(first_index, second_index)); 
   }

private:
   template<typename ITERATOR>
       std::size_t set_dimensions(ITERATOR begin, ITERATOR end)
   {
      // first calculate amount of memory needed in bytes
      const auto s = std::distance(begin, end);
      offsets_.clear();
      offsets_.reserve(s);
      offsets_.push_back(0);
      for(auto it=begin; it!=end; ++it) {
         assert(*it >= 0);
         offsets_.push_back( offsets_.back() + *it );
      }
      return offsets_.back();
   }

   template<typename ITERATOR>
   static std::vector<std::size_t> compute_offsets(ITERATOR begin, ITERATOR end)
   {
      // first calculate amount of memory needed in bytes
      const auto s = std::distance(begin, end);
      std::vector<std::size_t> offsets;
      offsets.reserve(s);
      offsets.push_back(0);
      for(auto it=begin; it!=end; ++it) {
         assert(*it >= 0);
         offsets.push_back( offsets.back() + *it );
      }
      return offsets;
   }

   std::vector<std::size_t> offsets_;
   std::vector<T> data_;
   // memory is laid out like this:
   // pointer[1], pointer[2], ... , pointer[dim1], pointer[end], data[1], data[2], ..., data[dim1] .
   // first pointers to the arrays in the second dimension are stored, then contiguously the data itself.
};

template<typename T>
void two_dim_variable_array<T>::clear()
{
    data_.clear();
    offsets_.clear();
}

template<typename T>
template<typename ITERATOR>
void two_dim_variable_array<T>::push_back(ITERATOR val_begin, ITERATOR val_end)
{
    offsets_.push_back(offsets_.back() + std::distance(val_begin, val_end));
    for(auto it=val_begin; it!=val_end; ++it)
    {
        data_.push_back(*it);
    }
}

} // namespace LPMP

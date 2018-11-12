#ifndef LPMP_SERIALIZE_HXX
#define LPMP_SERIALIZE_HXX

#include <array>
#include <vector>
#include "vector.hxx"
#include <bitset>
#include <cstring>

namespace LPMP {

  // to do: use perfect forwarding where applicable

template<typename T>
class binary_data {
public:
   binary_data(T* const _p, const INDEX _no_elements) : pointer(_p), no_elements(_no_elements) {}
   T* const pointer;
   const INDEX no_elements; 
};

class allocate_archive {
public:
  // for plain data
  template<typename T>
  static typename std::enable_if<std::is_arithmetic<T>::value,INDEX>::type
  serialize(const T&)
  {
    return sizeof(T); 
  }
  
  // for arrays
  template<typename T>
  static INDEX serialize(const T*, const INDEX s)
  {
    return sizeof(T)*s;
  }
  template<typename T>
  static INDEX serialize( const binary_data<T> b )
  {
     return serialize(b.pointer, b.no_elements); 
  }

  // for std::array<T>
  template<typename T, std::size_t N>
  static INDEX serialize(const std::array<T,N>&)
  {
    return sizeof(T)*N;
  }

  // for vector<T>
  template<typename T>
  static INDEX serialize(const vector<T>& v)
  {
     return serialize(v.begin(), v.size());
  }

  template<typename T>
  static INDEX serialize(const matrix<T>& m)
  {
     return m.size()*sizeof(T);
  }

  // for std::vector<T>
  template<typename T>
  static INDEX serialize(const std::vector<T>& v)
  {
     return serialize(v.data(), v.size());
  }

  // for std::bitset<N>
  template<std::size_t N>
  static INDEX serialize(const std::bitset<N>& v)
  {
     return sizeof(std::bitset<N>);
  }

  template<typename... T_REST>
  void operator()(T_REST&&... types)
  {}
  template<typename T, typename... T_REST>
  void operator()(T&& t, T_REST&&... types)
  {
     size_in_bytes_ += serialize(t);
     (*this)(types...);
  }

  INDEX size() const
  {
    return size_in_bytes_; 
  }

private:
  INDEX size_in_bytes_ = 0;
};

class serialization_archive {
public:
   serialization_archive() {}

  template<typename ITERATOR, typename SERIALIZATION_FUN>
  serialization_archive(ITERATOR begin, ITERATOR end, SERIALIZATION_FUN serialization_fun)
  {
    allocate_archive s;
    for(auto it=begin; it!=end; ++it) {
      serialization_fun(*it,s);
    }

    const INDEX size_in_bytes = s.size();

    archive_ = new char[size_in_bytes];
    end_ = archive_ + size_in_bytes;
    assert(archive_ != nullptr);
    cur_ = archive_;
  }

  serialization_archive(const allocate_archive& a)
  {
     const INDEX size_in_bytes = a.size();

     archive_ = new char[size_in_bytes];
     assert(archive_ != nullptr);
     end_ = archive_ + size_in_bytes;
     cur_ = archive_;
  }

  serialization_archive(const serialization_archive& o)
  {
     const INDEX size_in_bytes = o.size();

     archive_ = new char[size_in_bytes];
     assert(archive_ != nullptr);
     end_ = archive_ + size_in_bytes;
     cur_ = archive_ + (o.cur_ - o.archive_);

     std::memcpy(archive_, o.archive_, size_in_bytes); 
  }

  serialization_archive(serialization_archive&& o)
  {
     archive_ = o.archive_;
     end_ = o.end_;
     cur_ = o.cur_;

     o.archive_ = nullptr;
     o.end_ = nullptr;
     o.cur_ = nullptr;
  }

  serialization_archive(const void* mem, INDEX size_in_bytes)
  {
     assert(mem != nullptr);
     archive_ = (char*) mem;
     cur_ = archive_;
     end_ = archive_ + size_in_bytes;
  }

  ~serialization_archive()
  {
     free_memory();
  }

  void aquire_memory(const INDEX size_in_bytes)
  {
     free_memory();
     archive_ = new char[size_in_bytes];
     end_ = archive_ + size_in_bytes;
     assert(archive_ != nullptr);
     cur_ = archive_;
  }

  void release_memory()
  {
      archive_ = nullptr;
      cur_ = nullptr;
      end_ = nullptr; 
  }
  void free_memory()
  {
     if(archive_ != nullptr) {
        delete[] archive_;
        archive_ = nullptr;
        cur_ = nullptr;
        end_ = nullptr;
     }
  }

  bool operator==(const serialization_archive& o) const
  {
     if(size() != o.size()) { 
        return false;
     }
     return std::memcmp(archive_, o.archive_, size()) == 0;
  }

  INDEX size() const 
  {
     return end_ - archive_;
  }

  char* cur_address() const
  {
    return cur_;
  }

  void advance(const SIGNED_INDEX bytes) 
  {
    cur_ += bytes;
    assert(cur_ >= archive_);
    assert(cur_ <= end_);
  }

  void reset_cur()
  {
     cur_ = archive_;
  } 

  char* begin() { return archive_; }
  char* end() { return end_; }

private:
  char* archive_ = nullptr;
  char* end_ = nullptr;
  char* cur_;
};

// write data into archive
class save_archive {
public:
   save_archive(serialization_archive& a) 
      : ar(a) 
   {
      ar.reset_cur();
   }

  // for arrays
  template<typename T>
  void serialize(const T* p, const INDEX size)
  {
      INDEX* s = (INDEX*) ar.cur_address();
      std::memcpy((void*) s, (void*) p, size*sizeof(T));
      const INDEX size_in_bytes = sizeof(T)*size;
      ar.advance(size_in_bytes);
  }
  template<typename T>
  void serialize( const binary_data<T> b )
  {
     serialize(b.pointer, b.no_elements); 
  }

  // for std::array<T,N>
  template<typename T, std::size_t N>
  void serialize(const std::array<T,N>& v)
  {
    T* s = (T*) ar.cur_address();
    for(auto x : v) {
       *s = x;
       ++s; 
    }
    const INDEX size_in_bytes = sizeof(T)*N;
    ar.advance(size_in_bytes);
  } 
 
  // for vector<T>
  template<typename T>
  void serialize(const vector<T>& v)
  {
     serialize(v.begin(), v.size());
  }

  template<typename T>
  void serialize(const matrix<T>& m)
  {
     T* s = (T*) ar.cur_address();
     for(std::size_t i=0; i<m.dim1(); ++i) {
         for(std::size_t j=0; j<m.dim2(); ++j) {
             *s = m(i,j);
             ++s;
         }
     }
     //for(auto it=m.begin(); it!=m.end(); ++it) {
     //  *s = *it;
     //  ++s; 
     //}
     const INDEX size_in_bytes = sizeof(T)*m.size();
     ar.advance(size_in_bytes);
  }

  // for std::vector<T>
  template<typename T>
  void serialize(const std::vector<T>& v)
  {
     serialize(v.data(), v.size());
  }

  // for std::bitset<N>
  template<std::size_t N>
  static INDEX serialize(const std::bitset<N>& v)
  {
    assert(false);
    return 0;
  }

   // for plain data
   template<typename T>
   typename std::enable_if<std::is_arithmetic<T>::value>::type
   serialize(T& t)
   {
      T* s = (T*) ar.cur_address();
      *s = t;
      ar.advance(sizeof(T));
   }


   // save multiple entries
  template<typename... T_REST>
  void operator()(T_REST&&... types)
  {}
  template<typename T, typename... T_REST>
  void operator()(T&& t, T_REST&&... types)
  {
     serialize(t);
     (*this)(types...);
  }

private:
  serialization_archive& ar; 
};

// write data from archive into objects
class load_archive
{
public:
   load_archive(serialization_archive& a) 
      : ar(a) 
   {
      ar.reset_cur(); 
   }

   // for arrays
  template<typename T>
  void serialize(T* pointer, const INDEX size)
  {
      T* s = (T*) ar.cur_address();
      std::memcpy((void*) pointer, (void*) s, size*sizeof(T));
      const INDEX size_in_bytes = sizeof(T)*size;
      ar.advance(size_in_bytes);
  }
  template<typename T>
  void serialize( binary_data<T> b )
  {
     serialize(b.pointer, b.no_elements); 
  }

  // for std::array<T,N>
  template<typename T, std::size_t N>
  void serialize(std::array<T,N>& v)
  {
    T* s = (T*) ar.cur_address();
    for(auto& x : v) {
       x = *s;
       ++s; 
    }
    const INDEX size_in_bytes = sizeof(T)*N;
    ar.advance(size_in_bytes);
  } 
 
  // for vector<T>
  template<typename T>
  void serialize(vector<T>& v)
  {
     serialize(v.begin(), v.size());
  }

  template<typename T>
  void serialize(matrix<T>& m)
  {
     T* s = (T*) ar.cur_address();
     for(auto it=m.begin(); it!=m.end(); ++it) {
       *it = *s;
       ++s; 
     }
     const INDEX size_in_bytes = sizeof(T)*m.size();
     ar.advance(size_in_bytes);
  }

  // for std::vector<T>
  template<typename T>
  void serialize(std::vector<T>& v)
  {
     serialize(v.data(), v.size());
  } 

  // for std::bitset<N>
  template<std::size_t N>
  static INDEX serialize(const std::bitset<N>& v)
  {
    assert(false);
    return 0;
  }


   // for plain data
   template<typename T>
   typename std::enable_if<std::is_arithmetic<T>::value>::type
   serialize(T& t)
   {
      T* s = (T*) ar.cur_address();
      t = *s;
      ar.advance(sizeof(T));
   }

   // save multiple entries
  template<typename... T_REST>
  void operator()(T_REST&&... types)
  {}
  template<typename T, typename... T_REST>
  void operator()(T&& t, T_REST&&... types)
  {
     serialize(t);
     (*this)(types...);
  }

private:
  serialization_archive& ar;
};

// add numeric values stored in archive to variables
// do arithmetic on values stored in archive
enum class operation {addition, subtraction, multiplication, division, set_to_value};

template<operation OPERATION>
class arithmetic_archive {
public:
   arithmetic_archive(const REAL val) : val_(val) {}

   // for arrays
   template<typename T, typename OP>
     void serialize(T* pointer, const INDEX size, OP op)
     {
       static_assert(std::is_same<T,float>::value || std::is_same<T,double>::value,"");
       for(INDEX i=0; i<size; ++i) {
         pointer[i] = op(pointer[i], T(val_));
       }
     }
   template<typename T, typename OP>
     void serialize( binary_data<T> b, OP op)
     {
       serialize(b.pointer, b.no_elements, op); 
     }

   // for std::array<T,N>
   template<typename T, std::size_t N, typename OP>
     void serialize(std::array<T,N>& v, OP op)
     {
       static_assert(std::is_same<T,float>::value || std::is_same<T,double>::value,"");
       for(auto& x : v) {
          x = op(x, T(val_));
       }
     } 

   // for vector<T>
   template<typename T, typename OP>
     void serialize(vector<T>& v, OP op)
     {
       serialize(v.begin(), v.size(), op);
     }

   template<typename T, typename OP>
     void serialize(matrix<T>& m, OP op)
     {
       static_assert(std::is_same<T,float>::value || std::is_same<T,double>::value,"");
       for(auto it=m.begin(); it!=m.end(); ++it) {
         *it = op(*it, T(val_));
       }
     }

   // for std::vector<T>
   template<typename T, typename OP>
     void serialize(std::vector<T>& v, OP op)
     {
       serialize(v.data(), v.size(), op);
     } 

   // for plain data
   template<typename T, typename OP>
     typename std::enable_if<std::is_arithmetic<T>::value>::type
     serialize(T& t, OP op)
     {
       static_assert(std::is_same<T,float>::value || std::is_same<T,double>::value,"");
       t = op(t, T(val_));
     }

   // save multiple entries
   template<typename... T_REST>
     void operator()(T_REST&&... types)
     {}
   template<typename T, typename... T_REST>
     void operator()(T&& t, T_REST&&... types)
     {
        if(OPERATION == operation::addition) {
           auto op = [](auto a, auto b) { return a + b; };
           serialize(t, op);
        } else if(OPERATION == operation::subtraction) {
           assert(false);
        } else if(OPERATION == operation::multiplication) {
           auto op = [](auto a, auto b) { return a * b; };
           serialize(t, op);
        } else if(OPERATION == operation::division) {
           auto op = [](auto a, auto b) { return a / b; };
           serialize(t, op);
        } 
        else if(OPERATION == operation::set_to_value) {
           auto op = [](auto a, auto b) { return b; };
           serialize(t, op);
        } else {
           assert(false);
        }

       (*this)(types...);
     }

private:
   REAL val_;
};


class addition_archive {
public:
   addition_archive(serialization_archive& a, const REAL scaling) 
      : ar(a),
      scaling_(scaling) 
   {
      ar.reset_cur(); 
   }

   // for arrays
   template<typename T>
     void serialize(T* pointer, const INDEX size)
     {
       static_assert(std::is_same<T,float>::value || std::is_same<T,double>::value,"");
       T* val = (T*) ar.cur_address();
       for(INDEX i=0; i<size; ++i) {
         pointer[i] += T(scaling_) * val[i]; 
       }
       const INDEX size_in_bytes = sizeof(T)*size;
       ar.advance(size_in_bytes);
     }
   template<typename T>
     void serialize( binary_data<T> b )
     {
       serialize(b.pointer, b.no_elements); 
     }

   // for std::array<T,N>
   template<typename T, std::size_t N>
     void serialize(std::array<T,N>& v)
     {
       static_assert(std::is_same<T,float>::value || std::is_same<T,double>::value,"");
       T* s = (T*) ar.cur_address();
       for(auto& x : v) {
         x += T(scaling_) * (*s);
         ++s; 
       }
       const INDEX size_in_bytes = sizeof(T)*N;
       ar.advance(size_in_bytes);
     } 

   // for vector<T>
   template<typename T>
     void serialize(vector<T>& v)
     {
       serialize(v.begin(), v.size());
     }

   template<typename T>
     void serialize(matrix<T>& m)
     {
       static_assert(std::is_same<T,float>::value || std::is_same<T,double>::value,"");
       T* val = (T*) ar.cur_address();
       for(auto it=m.begin(); it!=m.end(); ++it) {
         *it += T(scaling_) * (*val); 
         ++val;
       }
       const INDEX size_in_bytes = sizeof(T)*m.size();
       ar.advance(size_in_bytes);
     }

   // for std::vector<T>
   template<typename T>
     void serialize(std::vector<T>& v)
     {
       serialize(v.data(), v.size());
     } 

  // for std::bitset<N>
  template<std::size_t N>
  static INDEX serialize(const std::bitset<N>& v)
  {
    assert(false);
    return 0;
  }


   // for plain data
   template<typename T>
     typename std::enable_if<std::is_arithmetic<T>::value>::type
     serialize(T& t)
     {
       static_assert(std::is_same<T,float>::value || std::is_same<T,double>::value,"");
       T* s = (T*) ar.cur_address();
       t += T(scaling_) * (*s);
       ar.advance(sizeof(T));
     }

   // save multiple entries
   template<typename... T_REST>
     void operator()(T_REST&&... types)
     {}
   template<typename T, typename... T_REST>
     void operator()(T&& t, T_REST&&... types)
     {
       serialize(t);
       (*this)(types...);
     }

private:
  serialization_archive& ar;
  const REAL scaling_;
};

} // end namespace LPMP
#endif // LPMP_SERIALIZE_HXX


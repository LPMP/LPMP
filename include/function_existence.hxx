#ifndef LPMP_FUNCTION_EXISTENCE_HXX
#define LPMP_FUNCTION_EXISTENCE_HXX

#include <type_traits>


// generates helper classes and constexpr function with which it is possible to detect existence and callability of member functions
// question: does this also help with inherited member functions?

#define LPMP_FUNCTION_EXISTENCE_CLASS(TESTER_NAME, MEMBER) \
template<typename, typename T> \
struct struct_##TESTER_NAME { \
   static_assert( \
         std::integral_constant<T, false>::value, \
         "Second template parameter needs to be of function type."); \
};  \
\
template<typename C, typename Ret, typename... Args> \
struct struct_##TESTER_NAME <C,Ret(Args...)> { \
private: \
   template<typename T> \
   static constexpr auto check(T*) \
   -> typename \
      std::is_same< \
      decltype( std::declval<T>() .MEMBER ( std::declval<Args>()...) ), \
      Ret \
         >::type { \
           return typename  std::is_same< decltype( std::declval<T>() .MEMBER ( std::declval<Args>()...) ), Ret >::type{}; } \
\
   template<typename> \
   static constexpr std::false_type check(...) { return std::false_type{}; } \
\
   typedef decltype(check<C>(nullptr)) type; \
\
public: \
   static constexpr bool value = type::value; \
}; \
template<class C, typename RET, typename... ARGS> \
constexpr static bool TESTER_NAME () { \
  return struct_##TESTER_NAME<C,RET(ARGS&...)>::value; \
} \




// generates helper class and constexpr function with which is is possible to detect existence and callability of member function that returns an object that one can write to
// do zrobienia: possibly also check for operator+= and operator-=
#define LPMP_ASSIGNMENT_FUNCTION_EXISTENCE_CLASS(TESTER_NAME, MEMBER) \
template<typename C, typename VALUE_TYPE, typename... Args> \
struct struct_##TESTER_NAME { \
private: \
   template<typename T> \
      static constexpr auto check(T*) \
      -> typename \
      std::is_same< \
      decltype( VALUE_TYPE (std::declval<T>(). MEMBER (std::declval<Args>()...) = std::declval<VALUE_TYPE>() ) ), \
      VALUE_TYPE \
         >::type { return typename std::is_same<bool,bool>::type{}; }; \
 \
   template<typename> \
      static constexpr std::false_type check(...) { return std::false_type{}; }; \
 \
   typedef decltype(check<C>(nullptr)) type; \
public: \
   static constexpr bool value = type::value; \
}; \
template<class C, typename VALUE_TYPE, typename... ARGS> \
constexpr static bool TESTER_NAME () { \
   return struct_##TESTER_NAME<C,VALUE_TYPE,ARGS&...>::value; \
} \


// generates helper class and constexpre function with which it is possible to detect existence of class existence
// do zrobienia: stupid name, but consistent with above. Possibly rename all of them
#define LPMP_CLASS_EXISTENCE_CLASS(TESTER_NAME,CLASS_NAME) \
template<typename C> \
struct struct_##TESTER_NAME { \
private: \
   template<typename T> \
      static constexpr auto check(T*) \
   -> typename std::is_class<typename T::CLASS_NAME>::type; \
   template<typename> \
      static constexpr std::false_type check(...); \
   typedef decltype(check<C>(nullptr)) type; \
public: \
   static constexpr bool value = type::value; \
}; \
template<typename C> \
constexpr static bool TESTER_NAME () { \
   return struct_##TESTER_NAME<C>::value; \
} \


#endif // LPMP_FUNCTION_EXISTENCE_HXX

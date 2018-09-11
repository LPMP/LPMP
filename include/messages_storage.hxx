#ifndef LPMP_MESSAGES_STORAGE_HXX
#define LPMP_MESSAGES_STORAGE_HXX

#include <vector>
#include <cassert>
#include "config.hxx"
#include "meta/meta.hpp"
#include "factor_container_interface.h"

namespace LPMP {

template<typename FMC>
class messages_storage {
public:
   struct message_info
   {
       FactorTypeAdapter* left;
       FactorTypeAdapter* right;
       message_passing_schedule::type mps;
   };

   message_info get_message(const std::size_t i) const { assert(i<messages_.size()); return messages_[i]; } 
   std::size_t number_of_messages() const { return messages_.size(); }

   template<typename MESSAGE_CONTAINER_TYPE, typename LEFT_FACTOR, typename RIGHT_FACTOR, typename... ARGS>
   MESSAGE_CONTAINER_TYPE* add_message(LEFT_FACTOR* l, RIGHT_FACTOR* r, ARGS... args);

   auto begin() { return messages_.begin(); }
   auto end() { return messages_.end(); }

private:
   template<typename MESSAGE_CONTAINER_TYPE>
   static constexpr std::size_t messages_tuple_index();

   std::vector<message_info> messages_; 

   struct vector_of_pointers {
      template<class T> 
         using invoke = typename std::vector<T*>;
   };

   using message_vector_list = meta::transform< typename FMC::MessageList, vector_of_pointers >;
   using message_storage_type = meta::apply<meta::quote<std::tuple>, message_vector_list>;

   message_storage_type messages_tuple_;
};

template<typename FMC>
template<typename MESSAGE_CONTAINER_TYPE, typename LEFT_FACTOR, typename RIGHT_FACTOR, typename... ARGS>
MESSAGE_CONTAINER_TYPE* messages_storage<FMC>::add_message(LEFT_FACTOR* l, RIGHT_FACTOR* r, ARGS... args)
{
   MESSAGE_CONTAINER_TYPE* m_l = l->template add_message<MESSAGE_CONTAINER_TYPE,Chirality::left>(r,std::forward<ARGS>(args)...);
   MESSAGE_CONTAINER_TYPE* m_r = r->template add_message<MESSAGE_CONTAINER_TYPE,Chirality::right>(l,std::forward<ARGS>(args)...);
   assert(m_l != nullptr || m_r != nullptr);

   l->set_left_msg(m_r);
   r->set_right_msg(m_l); 

   MESSAGE_CONTAINER_TYPE* m = (m_l != nullptr) ? m_l : m_r;
   assert(m != nullptr && (m == m_l || m == m_r));
   messages_.push_back({l,r, m->get_message_passing_schedule()}); // TODO get_message_passing_schedule should be accessible from MESSAGE_CONTAINER_TYPE

   constexpr auto msg_idx = messages_tuple_index<MESSAGE_CONTAINER_TYPE>();
   std::get<msg_idx>(messages_tuple_).push_back(m);
   return m;
}

template<typename FMC>
template<typename MESSAGE_CONTAINER_TYPE>
constexpr std::size_t messages_storage<FMC>::messages_tuple_index()
{
   constexpr auto n = meta::find_index<typename FMC::MessageList, MESSAGE_CONTAINER_TYPE>::value;
   static_assert(n < meta::size<typename FMC::MessageList>::value);
   return n;
}

} // namespace LPMP

#endif // LPMP_MESSAGES_STORAGE_HXX

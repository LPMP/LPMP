#pragma once

#include "pegtl.hh"

// elementary grammar rules for PEGTL

namespace LPMP {

namespace parsing {

struct mand_whitespace : pegtl::plus< pegtl::blank > {}; 
struct opt_whitespace : pegtl::star< pegtl::blank > {}; 
struct opt_invisible : pegtl::star< pegtl::sor< pegtl::blank, pegtl::eol > > {};
struct mand_invisible : pegtl::plus< pegtl::sor< pegtl::blank, pegtl::eol > > {};
struct positive_integer : pegtl::plus< pegtl::digit > {};

struct real_number_standard : pegtl::sor<
                              pegtl::seq< pegtl::opt< pegtl::one<'+','-'> >, pegtl::plus<pegtl::digit>, pegtl::opt< pegtl::seq< pegtl::string<'.'>, pegtl::star<pegtl::digit> > > >,
                              pegtl::string<'I','n','f'>,
                              pegtl::string<'i','n','f'>
                              > {}; 
struct real_number_smaller1 : pegtl::seq< pegtl::opt< pegtl::one<'+','-'> >, pegtl::string<'.'>, pegtl::plus< pegtl::digit > > {};
struct real_number_exponential : pegtl::seq< pegtl::opt< pegtl::one<'+','-'> >, pegtl::star< pegtl::digit >, pegtl::opt<pegtl::seq< pegtl::string<'.'>, pegtl::star< pegtl::digit>>>, pegtl::string<'e'>, pegtl::opt< pegtl::one<'+','-'> >, pegtl::plus< pegtl::digit > > {};
struct real_number : pegtl::sor<real_number_exponential, real_number_standard, real_number_smaller1> {};

struct vector : pegtl::seq< pegtl::string<'['>, opt_whitespace, real_number, pegtl::star< pegtl::seq< mand_whitespace, real_number > >, opt_whitespace, pegtl::string<']'> > {};

} // namespace parsing

} // namespace LPMP

// TODO: test the above definitions with unit testing whether they accept and reject what they are supposed to
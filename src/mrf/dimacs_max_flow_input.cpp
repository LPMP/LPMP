#include <fstream>
#include <cassert>
#include <stdexcept>
#include <limits>
#include "pegtl.hh"
#include "pegtl_parse_rules.h"
#include "mrf/dimacs_max_flow_input.h"

namespace LPMP {

namespace dimacs_max_flow_input {

// import basic parsers
using parsing::opt_whitespace;
using parsing::mand_whitespace;
using parsing::positive_integer;
using parsing::real_number;

struct empty_line : pegtl::seq<opt_whitespace, pegtl::eolf> {};
struct comment_line : pegtl::seq< opt_whitespace, pegtl::string<'c'>, pegtl::until<pegtl::eolf> > {};
struct ignore_line : pegtl::sor<comment_line, empty_line> {};

struct init_line : pegtl::seq< opt_whitespace, pegtl::string<'p'>, mand_whitespace, pegtl::string<'m','a','x'>, mand_whitespace, positive_integer, mand_whitespace, positive_integer, opt_whitespace, pegtl::eol> {};
struct source_node_line : pegtl::seq< opt_whitespace, pegtl::string<'n'>, mand_whitespace, positive_integer, mand_whitespace, pegtl::string<'s'>, opt_whitespace, pegtl::eol > {};
struct terminal_node_line : pegtl::seq< opt_whitespace, pegtl::string<'n'>, mand_whitespace, positive_integer, mand_whitespace, pegtl::string<'t'>, opt_whitespace, pegtl::eol > {};
struct arc_line : pegtl::seq< opt_whitespace, pegtl::string<'a'>, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number, opt_whitespace, pegtl::eolf > {};

struct grammar : pegtl::seq<
                 pegtl::star<ignore_line>,
                 init_line,
                 pegtl::star<ignore_line>,
                 source_node_line,
                 pegtl::star<ignore_line>,
                 terminal_node_line,
                 pegtl::star< pegtl::sor<ignore_line,arc_line> >
                 > {};

template< typename Rule >
struct action
: pegtl::nothing< Rule > {};

template<> struct action< init_line > {
    template<typename Input>
        static void apply(const Input& in, max_flow_instance& mc) 
        {
        }
};
template<> struct action< source_node_line > {
    template<typename Input>
        static void apply(const Input& in, max_flow_instance& mc) 
        {
            std::istringstream iss(in.string());
            unsigned char c; iss >> c;
            assert(c == 'n');
            iss >> mc.source;
            mc.source--;
            iss >> c;
            assert(c == 's');
        }
};
template<> struct action< terminal_node_line > {
    template<typename Input>
        static void apply(const Input& in, max_flow_instance& mc) 
        {
            std::istringstream iss(in.string());
            unsigned char c; iss >> c;
            assert(c == 'n');
            iss >> mc.terminal;
            mc.terminal--;
            iss >> c;
            assert(c == 't');
        }
};
template<> struct action< arc_line > {
    template<typename Input>
        static void apply(const Input& in, max_flow_instance& mc) 
        {
            std::istringstream iss(in.string());
            unsigned char a;
            iss >> a;
            assert(a == 'a');

            std::size_t i;
            iss >> i;
            assert(i >= 1);
            i--; 

            std::size_t j;
            iss >> j;
            assert(j >= 1);
            j--;

            double capacity = -std::numeric_limits<double>::infinity();
            iss >> capacity;
            assert(capacity >= 0.0); 

            mc.add_arc(i,j,capacity);
        }
};

max_flow_instance parse_file(const std::string& filename)
{
       max_flow_instance input;
       pegtl::file_parser problem(filename);
       const bool read_success = problem.parse< grammar, action >(input);
       if(!read_success) {
           throw std::runtime_error("could not read file " + filename);
       }

       assert(input.source == 0);
       assert(input.terminal == 1);

       return input;
}

max_flow_instance parse_string(const std::string& string)
{
       max_flow_instance input;
       const bool read_success = pegtl::parse<grammar, action>(string,"",input);
       if(!read_success) {
           throw std::runtime_error("could not read string");
       }

       assert(input.source == 0);
       assert(input.terminal == 1);

       return input;
}

} // namespace dimacs_max_flow_input

} // namespace LPMP

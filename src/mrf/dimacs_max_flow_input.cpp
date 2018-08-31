#include <fstream>
#include <cassert>
#include <stdexcept>
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

struct init_line : pegtl::seq< opt_whitespace, pegtl::string<'p'>, mand_whitespace, pegtl::string<'m','a','x'>, mand_whitespace, positive_integer, mand_whitespace, positive_integer, opt_whitespace, pegtl::eol> {};
struct source_node_line : pegtl::seq< opt_whitespace, pegtl::string<'n'>, mand_whitespace, positive_integer, mand_whitespace, pegtl::string<'s'>, opt_whitespace, pegtl::eol > {};
struct terminal_node_line : pegtl::seq< opt_whitespace, pegtl::string<'n'>, mand_whitespace, positive_integer, mand_whitespace, pegtl::string<'t'>, opt_whitespace, pegtl::eol > {};
struct arc_line : pegtl::seq< opt_whitespace, pegtl::string<'a'>, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number, opt_whitespace, pegtl::eolf > {};

struct grammar : pegtl::seq<
                 init_line,
                 source_node_line,
                 terminal_node_line,
                 //pegtl::opt<source_node_line, terminal_node_line>,
                 //pegtl::opt<source_node_line, terminal_node_line>,
                 pegtl::star<arc_line>
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
            typename max_flow_instance::capacitated_arc arc;
            unsigned char a;
            iss >> a;
            assert(a == 'a');

            iss >> arc[0];
            assert(arc[0] >= 1);
            arc[0]--;
            
            iss >> arc[1];
            assert(arc[1] >= 1);
            arc[1]--;
            
            iss >> arc.capacity;
            assert(arc.capacity >= 0.0); 

            mc.arcs.push_back(arc);
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

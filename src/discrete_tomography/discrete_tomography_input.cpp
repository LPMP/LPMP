#include "discrete_tomography/discrete_tomography_input.h"
#include "mrf/mrf_uai_input.h"
#include "pegtl.hh"
#include "pegtl_parse_rules.h"

namespace LPMP {

namespace discrete_tomography_UAI_input
{

using parsing::mand_whitespace;
using parsing::opt_whitespace;
using parsing::positive_integer;
using parsing::real_number;

struct ProjectionPreamble : pegtl::string<'P', 'R', 'O', 'J', 'E', 'C', 'T', 'I', 'O', 'N', 'S'> {};
struct ProjectionVector : pegtl::seq<pegtl::string<'('>, opt_whitespace, real_number, opt_whitespace, pegtl::star<pegtl::string<','>, opt_whitespace, real_number, opt_whitespace>, opt_whitespace, pegtl::string<')'>> {};
struct first_positive_integer : positive_integer {};
struct ProjectionLine : pegtl::seq<first_positive_integer, pegtl::plus<opt_whitespace, pegtl::string<'+'>, opt_whitespace, positive_integer>, opt_whitespace, pegtl::string<'='>, opt_whitespace, ProjectionVector> {};

// projection grammar
struct grammar : pegtl::seq<
                     pegtl::until<ProjectionPreamble>,
                     pegtl::star<pegtl::sor<pegtl::seq<opt_whitespace, pegtl::eol>, ProjectionLine>>,
                     pegtl::eof>
{};

template <typename Rule>
struct action : pegtl::nothing<Rule>
{};

template <>
struct action<pegtl::string<'('>>
{
    template <typename INPUT>
    static void apply(const INPUT &in, discrete_tomography_instance &instance)
    {
        instance.projection_costs.push_back({});
    }
};

template <>
struct action<real_number>
{
    template <typename INPUT>
    static void apply(const INPUT &in, discrete_tomography_instance &instance)
    {
        instance.projection_costs.back().push_back(std::stod(in.string()));
    }
};

template <>
struct action<first_positive_integer>
{
    template <typename INPUT>
    static void apply(const INPUT &in, discrete_tomography_instance &instance)
    {
        instance.projection_variables.push_back({});
        instance.projection_variables.back().push_back(std::stoul(in.string()));
    }
};

template <>
struct action<positive_integer>
{
    template <typename INPUT>
    static void apply(const INPUT &in, discrete_tomography_instance &instance)
    {
        instance.projection_variables.back().push_back(std::stoul(in.string()));
    }
};

discrete_tomography_instance parse_file(const std::string &filename)
{
    discrete_tomography_instance instance;
    instance.mrf = LPMP::mrf_uai_input::parse_file(filename);

    pegtl::file_parser problem(filename);

    const bool ret = problem.parse<grammar, action>(instance);
    if (ret != true)
        throw std::runtime_error("could not read projection constraints for discrete tomography");
    assert(instance.projection_variables.size() == instance.projection_costs.size());
    instance.propagate_projection_costs();
    return instance;
}

discrete_tomography_instance parse_string(const std::string &input)
{
    discrete_tomography_instance instance;
    throw std::runtime_error("not implemented yet");
    instance.propagate_projection_costs();
    return instance; 
}


} // namespace DiscreteTomographyTextInput

} // namespace LPMP
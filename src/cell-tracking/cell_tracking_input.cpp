#include "cell-tracking/cell_tracking_input.h"
#include "pegtl.hh"
#include "pegtl_parse_rules.h"
#include <iostream>

namespace LPMP {

namespace cell_tracking_parser_2d {

   // import basic parsers
   using parsing::opt_whitespace;
   using parsing::mand_whitespace;
   using parsing::opt_invisible;
   using parsing::mand_invisible;
   using parsing::positive_integer;
   using parsing::real_number;

   struct my_eolf : pegtl::sor<pegtl::eol, pegtl::eof> {}; // TODO: revert to usual pegtl::eolf
   struct comment_line : pegtl::seq< opt_whitespace, pegtl::sor< pegtl::seq< pegtl::string<'#'>, pegtl::until<pegtl::eolf> >, pegtl::eolf> > {};

   // H timestep hypothesis_id {detection cost} (x_coord, lower coord) 
   struct cell_detection_hypothesis : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number, opt_whitespace, pegtl::string<'('>, opt_whitespace, real_number, opt_whitespace, pegtl::string<','>, opt_whitespace, real_number, opt_whitespace, pegtl::string<')'> > {};
   struct cell_detection_hypothesis_line : pegtl::seq< pegtl::string<'H'>, opt_whitespace, cell_detection_hypothesis, opt_whitespace, pegtl::eolf > {};

   struct construct_cell_indices : pegtl::success {};
   
   // appearance cost: appearance_id hypothesis_id cost
   struct appearance : pegtl::seq<positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {};
   struct appearance_line : pegtl::seq< pegtl::string<'A','P','P'>, mand_whitespace, appearance, opt_whitespace, pegtl::eolf > {};
   // disappearance cost: disappearance_id hypothesis_id cost
   struct disappearance : pegtl::seq<positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {};
   struct disappearance_line : pegtl::seq< pegtl::string<'D','I','S','A','P','P'>, mand_whitespace, disappearance, opt_whitespace, my_eolf > {};
   // CONFSET hyp_nr1 hyp_nr2 ... <= 1
   struct conflict : pegtl::seq< positive_integer, pegtl::star< pegtl::seq< opt_whitespace, pegtl::string<'+'>, opt_whitespace, positive_integer> > > {};
   struct conflict_line : pegtl::seq< pegtl::string<'C','O','N','F','S','E','T'>, opt_whitespace, conflict, opt_whitespace, pegtl::string<'<','='>, opt_whitespace, pegtl::string<'1'>, opt_whitespace, my_eolf> {};
   // identifier, first hypothesis id, second hypothesis id, cost
   struct mapping : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {}; 
   struct mapping_line : pegtl::seq< pegtl::string<'M','O','V','E'>, opt_whitespace, mapping, opt_whitespace, my_eolf> {};
   // identifier, left hypothesis id, right hypothesid 1, right hypothesis id 2, cost
   struct division : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {}; 
   struct division_line : pegtl::seq< pegtl::string<'D','I','V'>, opt_whitespace, division, opt_whitespace, my_eolf> {};

   // first come cell detection hypotheses, then {appearances, disappearances, mappings and divisions}, folloed by conflicts
   struct grammar : pegtl::seq<
                        pegtl::star<pegtl::sor<
                            cell_detection_hypothesis_line,
                            comment_line>>,

                        construct_cell_indices,

                        pegtl::star<pegtl::sor<
                            appearance_line,
                            disappearance_line,
                            mapping_line,
                            division_line,
                            comment_line>>,

                        pegtl::star<pegtl::sor<
                            comment_line,
                            conflict_line>>,

                        pegtl::eof>
   {
   };

   template< typename Rule >
     struct action
     : pegtl::nothing< Rule > {};

   template<> struct action< cell_detection_hypothesis > {
     template<typename INPUT>
     static void apply(const INPUT & in, cell_tracking_instance& input)
     {
       //std::cout << "cell detection: " << in.string() << "\n";
       std::stringstream s(in.string());
       std::size_t timestep; s >> timestep;
       std::size_t hypothesis_id; s >> hypothesis_id; 
       double detection_cost; s >> detection_cost;
       char opening_bracket; s >> opening_bracket; assert(opening_bracket == '(');
       double x; s >> x;
       char comma; s >> comma; assert(comma == ',');
       double y; s >> y;
       char closing_bracket; s >> closing_bracket; assert(closing_bracket == ')');

       //typename cell_tracking_instance::cell_detection c;
       //c.timestep = timestep;
       //c.cell_number = hypothesis_id;
       //c.detection_cost = detection_cost;

//       input.cell_detections.push_back(c);
       input.add_cell_detection(timestep, hypothesis_id, detection_cost);
     }
   };

   template<> struct action< construct_cell_indices > {
       template<typename INPUT>
           static void apply(const INPUT & in, cell_tracking_instance& input)
           {
             //std::cout << "construct cell indices\n";
             auto &cell_detections = input.cell_detections;

             assert(cell_detections.size() > 0);
             for(std::size_t i=0; i<cell_detections.size(); ++i)
             {
               assert(cell_detections[i].cell_number == i || cell_detections[i].cell_number == std::numeric_limits<std::size_t>::max());
             }

             //const auto nr_timesteps = cell_detections.back().timestep;
             //input.timestep_to_first_cell_index.resize(nr_timesteps + 1);
             //input.timestep_to_first_cell_index[0] = 0;
             //std::size_t current_timestep = 0;

             //for (std::size_t i = 0; i < cell_detections.size(); ++i)
             //{
             //  if (cell_detections[i].timestep > current_timestep)
             //  {
             //    const std::size_t next_timestep = cell_detections[i].timestep;
             //    for (std::size_t t = current_timestep + 1; t <= next_timestep; ++t)
             //    {
             //      input.timestep_to_first_cell_index[t] = i;
             //    }
             //  } 
             //  } 
           }
   };

   template<> struct action< appearance > {
     template<typename INPUT>
       static void apply(const INPUT & in, cell_tracking_instance& input)
       {
         //std::cout << "appearance: " << in.string() << "\n";
         std::stringstream s(in.string());
         std::size_t timestep; s >> timestep;
         std::size_t hypothesis_id; s >> hypothesis_id;
         double cost; s >> cost;
         //const auto cell_index = input.cell_index(timestep, hypothesis_id);
         assert(input.cell_detections[hypothesis_id].cell_number == hypothesis_id);
         input.cell_detections[hypothesis_id].appearance_cost = cost;
       }
   };

   template<> struct action< disappearance > {
     template<typename INPUT>
       static void apply(const INPUT & in, cell_tracking_instance& input)
       {
         //std::cout << "disappearance: " << in.string() << "\n";
         std::stringstream s(in.string());
         std::size_t timestep; s >> timestep;
         std::size_t hypothesis_id; s >> hypothesis_id;
         double cost; s >> cost;
         //const auto cell_index = input.cell_index(timestep, hypothesis_id);
         assert(input.cell_detections[hypothesis_id].cell_number == hypothesis_id);
         input.cell_detections[hypothesis_id].disappearance_cost = cost;
       }
   };

   template<> struct action< conflict > {
     template<typename INPUT>
     static void apply(const INPUT & in, cell_tracking_instance& input)
     {
       input.conflict_element_bounds.push_back(input.conflict_cells.size());

       std::stringstream s(in.string()); 
       std::size_t hypothesis_id; s >> hypothesis_id; 
       //input.conflict_cells.push_back( input.cell_index(timestep, hypothesis_id) );
       input.conflict_cells.push_back( hypothesis_id );
       char plus;
       while(s >> plus) {
         assert(plus == '+');
         std::size_t hypothesis_id; s >> hypothesis_id; 
         //input.conflict_cells.push_back( input.cell_index(timestep, hypothesis_id) );
         input.conflict_cells.push_back( hypothesis_id );
       }
     }
   };

   template<> struct action< mapping > {
     template<typename INPUT>
     static void apply(const INPUT & in, cell_tracking_instance& input)
     {
       //std::cout << "move line = " << in.string() << "\n";
       std::stringstream s(in.string());
       std::size_t id; s >> id;
       std::size_t cell_prev; s >> cell_prev;
       std::size_t cell_next; s >> cell_next;
       double cost; s >> cost;

       const auto cell_outgoing = cell_prev; //input.cell_index(timestep_prev, cell_prev);
       const auto cell_incoming = cell_next; //input.cell_index(timestep_next, cell_next);
       input.cell_transitions.push_back({cell_outgoing, cell_incoming, cost});

       //input.cell_detections[cell_outgoing].nr_outgoing_transition_edges++;
       //input.cell_detections[cell_incoming].nr_incoming_transition_edges++;
     }
   };

   template<> struct action< division > {
     template<typename INPUT>
     static void apply(const INPUT & in, cell_tracking_instance& input)
     {
       std::stringstream s(in.string());
       std::size_t id; s >> id;
       std::size_t cell_prev; s >> cell_prev;
       std::size_t cell_next_1; s >> cell_next_1;
       std::size_t cell_next_2; s >> cell_next_2;
       double cost; s >> cost;

       const auto cell_outgoing = cell_prev; //input.cell_index(timestep_prev, cell_prev);
       const auto cell_incoming_1 = cell_next_1; //input.cell_index(timestep_next_1, cell_next_1);
       const auto cell_incoming_2 = cell_next_2; //input.cell_index(timestep_next_2, cell_next_2);
       input.cell_divisions.push_back({cell_outgoing, cell_incoming_1, cell_incoming_2, cost});

       //input.cell_detections[cell_outgoing].nr_outgoing_division_edges++;
       //input.cell_detections[cell_incoming_1].nr_incoming_division_edges++;
       //input.cell_detections[cell_incoming_2].nr_incoming_division_edges++;
     }
   };


   cell_tracking_instance parse_file(const std::string& filename)
   {
       cell_tracking_instance input;
       pegtl::file_parser problem(filename);
       const bool success = problem.parse< grammar, action >(input); 
       if(!success)
         throw std::runtime_error("Could not read cell tracking input file " + filename);
       return input;
   }

   cell_tracking_instance parse_string(const std::string &filename)
   {
     cell_tracking_instance input;
     throw std::runtime_error("not implemented yet");
     return input;
   }

} // end namespace cell_tracking_parser_2d

} // namespace LPMP
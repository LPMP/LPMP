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
   struct comment_line : pegtl::seq< opt_whitespace, pegtl::sor< pegtl::seq< pegtl::string<'#'>, pegtl::until<my_eolf> >, my_eolf> > {};

   struct division_distance : pegtl::seq< positive_integer> {};
   struct division_distance_line : pegtl::seq< opt_whitespace, pegtl::string<'D','I','V','I','S','I','O','N',' ','D','I','S','T','A','N','C','E'>, opt_whitespace, pegtl::string<'='>, opt_whitespace, division_distance, opt_whitespace, my_eolf > {};
   
   // H timestep hypothesis_id {detection cost} (x_coord, lower coord) 
   struct cell_detection_hypothesis : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number, opt_whitespace, pegtl::string<'('>, opt_whitespace, real_number, opt_whitespace, pegtl::string<','>, opt_whitespace, real_number, opt_whitespace, pegtl::string<')'> > {};
   struct cell_detection_hypothesis_line : pegtl::seq< pegtl::string<'H'>, opt_whitespace, cell_detection_hypothesis, opt_whitespace, my_eolf > {};
   struct construct_cell_indices : pegtl::seq< pegtl::any > {};
   
   // appearance cost: timestep hypothesis_id cost
   struct appearance : pegtl::seq<positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {};
   struct appearance_line : pegtl::seq< pegtl::string<'A','P','P'>, mand_whitespace, appearance, opt_whitespace, my_eolf > {};
   // disappearance cost: timestep hypothesis_id cost
   struct disappearance : pegtl::seq<positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {};
   struct disappearance_line : pegtl::seq< pegtl::string<'D','I','S','A','P','P'>, mand_whitespace, disappearance, opt_whitespace, my_eolf > {};
   // CONFSET timestep_1 hyp_no1 + timestep hyp_no2 ... <= 1
   struct conflict : pegtl::seq< positive_integer, mand_whitespace, positive_integer, pegtl::star< pegtl::seq< opt_whitespace, pegtl::string<'+'>, opt_whitespace, positive_integer, mand_whitespace, positive_integer> > > {};
   struct conflict_line : pegtl::seq< pegtl::string<'C','O','N','F','S','E','T'>, opt_whitespace, conflict, opt_whitespace, pegtl::string<'<','='>, opt_whitespace, pegtl::string<'1'>, opt_whitespace, my_eolf> {};
   // timestep 1, first hypothesis id, timestep 1, second hypothesis id, cost
   struct mapping : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {}; 
   struct mapping_line : pegtl::seq< pegtl::string<'M','O','V','E'>, opt_whitespace, mapping, opt_whitespace, my_eolf> {};
   // timestep 1, left hypothesis id, right timestep 1, right hypothesid 1, right timestep 2, right hypothesis id 2, cost
   struct division : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {}; 
   struct division_line : pegtl::seq< pegtl::string<'D','I','V'>, opt_whitespace, division, opt_whitespace, my_eolf> {};


   // first come cell detection hypotheses, then {appearances, disappearances, mappings and divisions}, folloed by conflicts
   struct grammar : pegtl::seq< 
                    pegtl::star< pegtl::sor<
                        comment_line,
                        cell_detection_hypothesis_line
                    > >,

                    construct_cell_indices,

                    pegtl::star< pegtl::sor<
                        comment_line,
                        appearance_line,
                        disappearance_line,
                        mapping_line,
                        division_line
                    > >,

                    pegtl::star< pegtl::sor<
                        comment_line,
                        conflict_line
                    > >
                    > {};

   struct grammar_division_distance : 
                    pegtl::seq< 
                    pegtl::star< comment_line >,
                    division_distance_line,
                    
                    pegtl::star<pegtl::sor<
                    comment_line,
                    cell_detection_hypothesis_line,
                    appearance_line,
                    disappearance_line,
                    conflict_line,
                    mapping_line,
                    division_line
                    >>> {};


   template< typename Rule >
     struct action
     : pegtl::nothing< Rule > {};

   template<> struct action< division_distance > {
     template<typename INPUT>
     static void apply(const INPUT & in, cell_tracking_input& i)
     {
       assert(false);
       //assert(i.division_distance == 1);
       //std::stringstream s(in.string()); 
       //s >> i.division_distance; 
       //assert(i.division_distance > 0);
     }
   };

   template<> struct action< cell_detection_hypothesis > {
     template<typename INPUT>
     static void apply(const INPUT & in, cell_tracking_input& input)
     {
       std::stringstream s(in.string());
       std::size_t timestep; s >> timestep;
       std::size_t hypothesis_id; s >> hypothesis_id; 
       double detection_cost; s >> detection_cost;
       char opening_bracket; s >> opening_bracket; assert(opening_bracket == '(');
       double x; s >> x;
       char comma; s >> comma; assert(comma == ',');
       double y; s >> y;
       char closing_bracket; s >> closing_bracket; assert(closing_bracket == ')');

       typename cell_tracking_input::cell_detection c;
       c.timestep = timestep;
       c.cell_number = hypothesis_id;
       c.detection_cost = detection_cost;

       input.cell_detections.push_back(c);
     }
   };

   template<> struct action< construct_cell_indices > {
       template<typename INPUT>
           static void apply(const INPUT & in, cell_tracking_input& input)
           {
               auto& cell_detections = input.cell_detections;
               
               assert(cell_detections.size() > 0);
               assert(std::is_sorted(cell_detections.begin(), cell_detections.end()));

               const auto no_timesteps = cell_detections.back().timestep;
               input.timestep_to_first_cell_index.resize(no_timesteps+1);
               input.timestep_to_first_cell_index[0] = 0;
               std::size_t current_timestep = 0;

               for(std::size_t i=0; i<cell_detections.size(); ++i) {
                   if(cell_detections[i].timestep > current_timestep) {
                       const std::size_t next_timestep = cell_detections[i].timestep;
                       for(std::size_t t=current_timestep+1; t<=next_timestep; ++t) {
                           input.timestep_to_first_cell_index[t] = i;
                       } 
                   } 
               } 
           }
   };

   template<> struct action< appearance > {
     template<typename INPUT>
       static void apply(const INPUT & in, cell_tracking_input& input)
       {
         std::stringstream s(in.string());
         std::size_t timestep; s >> timestep;
         std::size_t hypothesis_id; s >> hypothesis_id;
         double cost; s >> cost;
         const auto cell_index = input.cell_index(timestep, hypothesis_id);
         input.cell_detections[cell_index].appearance_cost = 0;
       }
   };

   template<> struct action< disappearance > {
     template<typename INPUT>
       static void apply(const INPUT & in, cell_tracking_input& input)
       {
         std::stringstream s(in.string());
         std::size_t timestep; s >> timestep;
         std::size_t hypothesis_id; s >> hypothesis_id;
         double cost; s >> cost;
         const auto cell_index = input.cell_index(timestep, hypothesis_id);
         input.cell_detections[cell_index].disappearance_cost = 0;
       }
   };

   template<> struct action< conflict > {
     template<typename INPUT>
     static void apply(const INPUT & in, cell_tracking_input& input)
     {
       input.conflict_element_bounds.push_back(input.conflict_cells.size());

       std::stringstream s(in.string()); 
       std::size_t timestep; s >> timestep; 
       std::size_t hypothesis_id; s >> hypothesis_id; 
       input.conflict_cells.push_back( input.cell_index(timestep, hypothesis_id) );
       char plus;
       while(s >> plus) {
         assert(plus == '+');
         std::size_t timestep; s >> timestep; 
         std::size_t hypothesis_id; s >> hypothesis_id; 
         input.conflict_cells.push_back( input.cell_index(timestep, hypothesis_id) );
       }
     }
   };

   template<> struct action< mapping > {
     template<typename INPUT>
     static void apply(const INPUT & in, cell_tracking_input& input)
     {
       std::stringstream s(in.string());
       std::size_t timestep_prev; s >> timestep_prev;
       std::size_t cell_prev; s >> cell_prev;
       std::size_t timestep_next; s >> timestep_next;
       std::size_t cell_next; s >> cell_next;
       double cost; s >> cost;

       const auto cell_outgoing = input.cell_index(timestep_prev, cell_prev);
       const auto cell_incoming = input.cell_index(timestep_next, cell_next);
       input.cell_transitions.push_back({cell_outgoing, cell_incoming, cost});

       input.cell_detections[cell_outgoing].no_outgoing_transition_edges++;
       input.cell_detections[cell_incoming].no_incoming_transition_edges++;
     }
   };

   template<> struct action< division > {
     template<typename INPUT>
     static void apply(const INPUT & in, cell_tracking_input& input)
     {
       std::stringstream s(in.string());
       std::size_t timestep_prev; s >> timestep_prev;
       std::size_t cell_prev; s >> cell_prev;
       std::size_t timestep_next_1; s >> timestep_next_1;
       std::size_t cell_next_1; s >> cell_next_1;
       std::size_t timestep_next_2; s >> timestep_next_2;
       std::size_t cell_next_2; s >> cell_next_2;
       double cost; s >> cost;

       const auto cell_outgoing = input.cell_index(timestep_prev, cell_prev);
       const auto cell_incoming_1 = input.cell_index(timestep_next_1, cell_next_1);
       const auto cell_incoming_2 = input.cell_index(timestep_next_2, cell_next_2);
       input.cell_divisions.push_back({cell_outgoing, cell_incoming_1, cell_incoming_2, cost});

       input.cell_detections[cell_outgoing].no_outgoing_division_edges++;
       input.cell_detections[cell_incoming_1].no_incoming_division_edges++;
       input.cell_detections[cell_incoming_2].no_incoming_division_edges++;
     }
   };


/*
   template<typename GRAMMAR>
   bool read_input(const std::string& filename, cell_tracking_input& i)
   {
       std::cout << "parsing " << filename << "\n";
       pegtl::file_parser problem(filename);
       const bool success = problem.parse< GRAMMAR, action >(i); 
       return success;
   }

   template<typename SOLVER>
   void construct_tracking_problem(input& i, SOLVER& s)
   {
      auto& cell_tracking_constructor = s.template GetProblemConstructor<0>();
      auto& lp = s.GetLP();

      std::transform(i.cell_detection_stat.begin(), i.cell_detection_stat.end(),
              std::back_inserter(cell_tracking_constructor.cumulative_sum_cell_detection_factors),
              [](auto timestep_array) { return timestep_array.size(); }
              );
      std::partial_sum(cell_tracking_constructor.cumulative_sum_cell_detection_factors.begin(), cell_tracking_constructor.cumulative_sum_cell_detection_factors.end(), cell_tracking_constructor.cumulative_sum_cell_detection_factors.begin());
      cell_tracking_constructor.cumulative_sum_cell_detection_factors.back() = 0;
      std::rotate(cell_tracking_constructor.cumulative_sum_cell_detection_factors.begin(), cell_tracking_constructor.cumulative_sum_cell_detection_factors.end()-1, cell_tracking_constructor.cumulative_sum_cell_detection_factors.end());
      assert(cell_tracking_constructor.cumulative_sum_cell_detection_factors[0] == 0);

      const auto no_detection_hypotheses = std::accumulate(i.cell_detection_stat.begin(), i.cell_detection_stat.end(), 0, [](auto sum, auto timestep_array) { return sum + timestep_array.size(); });
      std::size_t no_cell_transitions = 0;
      std::size_t no_cell_divisions = 0;
      for(const auto& timestep_array : i.cell_detection_stat) {
          for(const auto& stat : timestep_array) {
              no_cell_transitions += stat.no_outgoing_transition_edges;
              no_cell_divisions += stat.no_outgoing_division_edges;
          }
      }
      cell_tracking_constructor.begin(lp, no_detection_hypotheses, no_cell_transitions, no_cell_divisions);

      //std::cout << "exclusion constraints disabled\n";
      cell_tracking_constructor.set_number_of_timesteps( i.cell_detection_stat.size() );
      for(std::size_t t=0; t<i.cell_detection_stat.size(); ++t) {
        for(std::size_t n=0; n<i.cell_detection_stat[t].size(); ++n) {
          const double detection_cost = i.cell_detection_stat[t][n].detection_cost;
          const double appearance_cost = i.cell_detection_stat[t][n].appearance_cost;
          const double disappearance_cost = i.cell_detection_stat[t][n].disappearance_cost;
          const std::size_t no_incoming_transition_edges = i.cell_detection_stat[t][n].no_incoming_transition_edges;
          const std::size_t no_outgoing_transition_edges = i.cell_detection_stat[t][n].no_outgoing_transition_edges;
          const std::size_t no_incoming_division_edges = i.cell_detection_stat[t][n].no_incoming_division_edges;
          const std::size_t no_outgoing_division_edges = i.cell_detection_stat[t][n].no_outgoing_division_edges;

          cell_tracking_constructor.add_detection_hypothesis( lp, t, n, detection_cost, appearance_cost, disappearance_cost, no_incoming_transition_edges, no_incoming_division_edges, no_outgoing_transition_edges, no_outgoing_division_edges);
        }
      }

      for(const auto& conflict_set : i.conflicts) {
        cell_tracking_constructor.add_exclusion_constraint(lp, conflict_set.begin(), conflict_set.end()); 
        //cell_tracking_constructor.register_exclusion_constraint(conflict_set.begin(), conflict_set.end() ); 
      }

      //auto tc = cell_tracking_constructor.init_transition_counter();
      // the order is important! Possibly change and treat transition and division edges separately
      for(auto& t : i.mappings) {
        const std::size_t timestep_prev = std::get<0>(t);
        const std::size_t prev_cell = std::get<1>(t);
        const std::size_t timestep_next = std::get<2>(t);
        const std::size_t next_cell = std::get<3>(t);
        const double cost = std::get<4>(t);
        cell_tracking_constructor.add_cell_transition( lp, timestep_prev, prev_cell, timestep_next, next_cell, cost);
      }
      for(auto& t : i.divisions) {
        const std::size_t timestep_prev = std::get<0>(t);
        const std::size_t prev_cell = std::get<1>(t);
        const std::size_t timestep_next_1 = std::get<2>(t);
        const std::size_t next_cell_1 = std::get<3>(t);
        const std::size_t timestep_next_2 = std::get<4>(t);
        const std::size_t next_cell_2 = std::get<5>(t);
        const double cost = std::get<6>(t);
        cell_tracking_constructor.add_cell_division( lp, timestep_prev, prev_cell, timestep_next_1, next_cell_1, timestep_next_2, next_cell_2, cost);
      }
      //for(std::size_t t=0; t<i.cell_detection_stat.size(); ++t) {
      //  for(std::size_t j=0; j<i.cell_detection_stat[t].size(); ++j) {
      //    assert( i.cell_detection_stat[t][j].no_outgoing_transition_edges + i.cell_detection_stat[t][j].no_outgoing_division_edges == tc.current_transition_no[t][j][1] + tc.current_division_no[t][j][1] );
      //    assert( i.cell_detection_stat[t][j].no_incoming_transition_edges + i.cell_detection_stat[t][j].no_incoming_division_edges == tc.current_transition_no[t][j][0] + tc.current_division_no[t][j][0] );
      //  }
      //}

      cell_tracking_constructor.end(lp);
   }
   */

   cell_tracking_input parse_file(const std::string& filename)
   {
       cell_tracking_input input;
       pegtl::file_parser problem(filename);
       const bool success = problem.parse< grammar, action >(input); 
       assert(success);
       return input;
   }
   cell_tracking_input parse_string(const std::string& filename)
   {}

} // end namespace cell_tracking_parser_2d

/*
namespace cell_tracking_parser_mother_machine {

  struct mother_machine_cell_tracking_input : cell_tracking_input {
    std::vector<std::vector<std::array<REAL,2>>> cell_position; // lower, upper
  };

   // import basic parsers
   using Parsing::opt_whitespace;
   using Parsing::mand_whitespace;
   using Parsing::opt_invisible;
   using Parsing::mand_invisible;
   using Parsing::positive_integer;
   using Parsing::real_number;

   struct my_eolf : pegtl::sor<pegtl::eol, pegtl::eof> {}; // do zrobienia: revert to usual pegtl::eolf
   struct comment_line : pegtl::seq< opt_whitespace, pegtl::sor< pegtl::seq< pegtl::string<'#'>, pegtl::until<my_eolf> >, my_eolf> > {};
   
   // t = number
   struct timestep : pegtl::seq< positive_integer > {};
   struct timestep_line : pegtl::seq< pegtl::string<'t'>, opt_whitespace, pegtl::string<'='>, opt_whitespace, timestep, opt_whitespace, my_eolf > {};

   // H hypothesis_number hypothesis_id {detection cost} {exit cost} (upper pos, lower pos) // upper < lower, inverted!
   struct cell_detection_hypothesis : pegtl::seq< positive_integer, opt_whitespace, positive_integer, opt_whitespace, real_number, opt_whitespace, real_number, opt_whitespace, pegtl::string<'('>, opt_whitespace, positive_integer, opt_whitespace, pegtl::string<','>, opt_whitespace, positive_integer, opt_whitespace, pegtl::string<')'> > {};
   struct cell_detection_hypothesis_line : pegtl::seq< pegtl::string<'H'>, opt_whitespace, cell_detection_hypothesis, opt_whitespace, my_eolf > {};
   
   // EC hyp_no1 + hyp_no2 ... <= 1
   struct exclusion : pegtl::seq< positive_integer, pegtl::star< pegtl::seq< opt_whitespace, pegtl::string<'+'>, opt_whitespace, positive_integer> > > {};
   struct exclusion_line : pegtl::seq< pegtl::string<'E','C'>, opt_whitespace, exclusion, opt_whitespace, pegtl::string<'<','='>, opt_whitespace, pegtl::string<'1'>, opt_whitespace, my_eolf> {};
   // timestep 1, first hypothesis id, timestep 1, second hypothesis id, cost
   struct mapping : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {}; 
   struct mapping_line : pegtl::seq< pegtl::string<'M','A'>, opt_whitespace, mapping, opt_whitespace, my_eolf> {};
   // timestep 1, left hypothesis id, timestep 2, two right hypothesis ids, cost
   struct division : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {}; 
   struct division_line : pegtl::seq< pegtl::string<'D','A'>, opt_whitespace, division, opt_whitespace, my_eolf> {};

   struct grammar : pegtl::seq< pegtl::star<pegtl::sor<
                    comment_line,
                    timestep_line, 
                    cell_detection_hypothesis_line,
                    exclusion_line,
                    mapping_line,
                    division_line
                    >>> {};

   template< typename Rule >
     struct action
     : pegtl::nothing< Rule > {};

   template<> struct action< pegtl::eof > {
     template<typename INPUT>
     static void apply(const INPUT & in, mother_machine_cell_tracking_input& i)
     {
       //std::cout << "kwaskwas\n";
       i.eof = true;
     }
   };

   template<> struct action< timestep > {
     template<typename INPUT>
     static void apply(const INPUT & in, mother_machine_cell_tracking_input& i)
     {
       const INDEX timestep = std::stoul(in.string());
       assert(i.cell_detection_stat.size() == timestep);
       i.cell_detection_stat.push_back({});
       i.cell_position.push_back({});
     }
   };

   template<> struct action< cell_detection_hypothesis > {
     template<typename INPUT>
     static void apply(const INPUT & in, mother_machine_cell_tracking_input& i)
     {
       std::stringstream s(in.string());
       INDEX number; s >> number;
       assert(number == i.cell_detection_stat.back().size());
       INDEX hypothesis_id; s >> hypothesis_id; // this number is not used
       REAL detection_cost; s >> detection_cost;
       REAL exit_cost; s >> exit_cost;
       char opening_bracket; s >> opening_bracket; assert(opening_bracket == '(');
       REAL upper_boundary; s >> upper_boundary;
       char comma; s >> comma; assert(comma == ',');
       REAL lower_boundary; s >> lower_boundary;
       char closing_bracket; s >> closing_bracket; assert(closing_bracket == ')');
       assert(upper_boundary < lower_boundary );

       i.cell_detection_stat.back().push_back({detection_cost,exit_cost,0,0,0,0});//, upper_boundary, lower_boundary));
       i.cell_position.back().push_back({lower_boundary, upper_boundary});

       //std::cout << "H: " << number << ", " << detection_cost << std::endl;
     }
   };

   template<> struct action< exclusion > {
     template<typename INPUT>
     static void apply(const INPUT & in, mother_machine_cell_tracking_input& i)
     {
       std::stringstream s(in.string());
       std::vector<std::array<INDEX,2>> idx;
       assert(i.cell_detection_stat.size() > 0);
       const INDEX timestep = i.cell_detection_stat.size()-1;
       INDEX number; s >> number; idx.push_back({timestep,number});
       //std::cout << "E: " << number << ", ";
       char plus;
       while(s >> plus) {
         assert(plus == '+');
         INDEX number; s >> number; idx.push_back({timestep,number}); 
         //std::cout << number << ", ";
       }
       i.conflicts.push_back(idx);
       //std::cout << " <= 1; " << i.exclusions_.size() << "\n";
     }
   };

   template<> struct action< mapping > {
     template<typename INPUT>
     static void apply(const INPUT & in, mother_machine_cell_tracking_input& i)
     {
       std::stringstream s(in.string());
       INDEX timestep_prev; s >> timestep_prev;
       INDEX cell_prev; s >> cell_prev;
       INDEX timestep_next; s >> timestep_next;
       INDEX cell_next; s >> cell_next;
       REAL cost; s >> cost;

       assert(timestep_prev+1 == timestep_next);
       assert(timestep_next < i.cell_detection_stat.size());
       assert(cell_prev < i.cell_detection_stat[timestep_prev].size());
       assert(cell_next < i.cell_detection_stat[timestep_next].size());

       i.cell_detection_stat[timestep_prev][cell_prev].no_outgoing_transition_edges++;
       i.cell_detection_stat[timestep_next][cell_next].no_incoming_transition_edges++;

       i.mappings.push_back( std::make_tuple(timestep_prev, cell_prev, timestep_next, cell_next, cost));

       //std::cout << "mapping: t = " << timestep_prev << ", h1 = " << cell_prev << ", h2 = " << cell_next << ", cost = " << cost << "\n";
     }
   };

   template<> struct action< division > {
     template<typename INPUT>
     static void apply(const INPUT & in, mother_machine_cell_tracking_input& i)
     {
       std::stringstream s(in.string());
       INDEX timestep_prev; s >> timestep_prev;
       INDEX cell_prev; s >> cell_prev;
       INDEX timestep_next; s >> timestep_next;
       INDEX cell_next_1; s >> cell_next_1;
       INDEX cell_next_2; s >> cell_next_2;
       REAL cost; s >> cost;

       assert(timestep_prev+1 == timestep_next);
       assert(cell_next_1 != cell_next_2);
       assert(timestep_next < i.cell_detection_stat.size());
       assert(cell_prev < i.cell_detection_stat[timestep_prev].size());
       assert(cell_next_1 < i.cell_detection_stat[timestep_next].size());
       assert(cell_next_2 < i.cell_detection_stat[timestep_next].size());

       i.cell_detection_stat[timestep_prev][cell_prev].no_outgoing_division_edges++;
       i.cell_detection_stat[timestep_next][cell_next_1].no_incoming_division_edges++;
       i.cell_detection_stat[timestep_next][cell_next_2].no_incoming_division_edges++;

       i.divisions.push_back( std::make_tuple(timestep_prev, cell_prev, timestep_next, cell_next_1, timestep_next, cell_next_2, cost));
       //std::cout << "division: t = " << timestep_prev << ", h1 = " << cell_prev << ", h1,1 = " << cell_next_1 << ", h2,2 = " << cell_next_2 << ", cost = " << cost << "\n";
     }
   };


   bool read_input(const std::string& filename, mother_machine_cell_tracking_input& i)
   {
      if(verbosity >= 1) {
        std::cout << "parsing " << filename << "\n";
      }
      pegtl::file_parser problem(filename);
      const bool success = problem.parse< grammar, action >(i); 
      //assert(i.eof);
      return success;
   }

   template<typename SOLVER>
   void construct_tracking_problem(mother_machine_cell_tracking_input& i, SOLVER& s)
   {
      auto& cell_tracking_constructor = s.template GetProblemConstructor<0>();
      auto& lp = s.GetLP();

      //begin(lp);

      //std::cout << "exclusion constraints disabled\n";
      cell_tracking_constructor.set_number_of_timesteps( i.cell_detection_stat.size() );
      for(INDEX t=0; t<i.cell_detection_stat.size(); ++t) {
        for(INDEX n=0; n<i.cell_detection_stat[t].size(); ++n) {
          const REAL detection_cost = i.cell_detection_stat[t][n].detection_cost;
          const REAL appearance_cost = i.cell_detection_stat[t][n].appearance_cost;
          const REAL disappearance_cost = i.cell_detection_stat[t][n].disappearance_cost;
          const INDEX no_incoming_transition_edges = i.cell_detection_stat[t][n].no_incoming_transition_edges;
          const INDEX no_outgoing_transition_edges = i.cell_detection_stat[t][n].no_outgoing_transition_edges;
          const INDEX no_incoming_division_edges = i.cell_detection_stat[t][n].no_incoming_division_edges;
          const INDEX no_outgoing_division_edges = i.cell_detection_stat[t][n].no_outgoing_division_edges;

          //assert(t != i.cell_detection_stat.size()-1 || exit_cost == 0.0);
          cell_tracking_constructor.add_detection_hypothesis( 
              lp, t, n, 
              detection_cost, appearance_cost, disappearance_cost,
              no_incoming_transition_edges, no_incoming_division_edges,
              no_outgoing_transition_edges, no_outgoing_division_edges 
              );
        }
        // check whether all cell hypotheses are covered by exclusion constraints (must be valid for mother machine)
        //std::vector<INDEX> all_indices;
        //for(const auto& detections : i.exclusions_[t]) {
        //  all_indices.insert(all_indices.end(), detections.begin(), detections.end());
        //}
        //std::sort(all_indices.begin(), all_indices.end());
        //all_indices.erase( std::unique( all_indices.begin(), all_indices.end()), all_indices.end() );
        //assert(all_indices.size() == i.cell_detection_stat[t].size());
        //
      }

      for(const auto& detections : i.conflicts) {
        cell_tracking_constructor.add_exclusion_constraint(lp, detections.begin(), detections.end() ); 
      }

      //auto tc = cell_tracking_constructor.init_transition_counter();
      for(auto& t : i.mappings) {
        const INDEX timestep_prev = std::get<0>(t);
        const INDEX prev_cell = std::get<1>(t);
        const INDEX timestep_next = std::get<2>(t);
        const INDEX next_cell = std::get<3>(t);
        const REAL cost = std::get<4>(t);
        cell_tracking_constructor.add_cell_transition( lp, timestep_prev, prev_cell, timestep_next, next_cell, cost);
      }
      for(auto& t : i.divisions) {
        const INDEX timestep_prev = std::get<0>(t);
        const INDEX prev_cell = std::get<1>(t);
        const INDEX timestep_next_1 = std::get<2>(t);
        const INDEX next_cell_1 = std::get<3>(t);
        const INDEX timestep_next_2 = std::get<4>(t);
        const INDEX next_cell_2 = std::get<5>(t);
        const REAL cost = std::get<6>(t);
        cell_tracking_constructor.add_cell_division( lp, timestep_prev, prev_cell, timestep_next_1, next_cell_1, timestep_next_2, next_cell_2, cost);
      }
      //for(INDEX t=0; t<i.cell_detection_stat.size(); ++t) {
      //  for(INDEX j=0; j<i.cell_detection_stat[t].size(); ++j) {
      //    assert( std::get<2>(i.cell_detection_stat[t][j]) == tc[t][j][0] );
      //    assert( std::get<3>(i.cell_detection_stat[t][j]) == tc[t][j][1] );
      //  }
      //}

       end(lp);
   }

   // we assume problem with single constructor
   template<typename SOLVER>
   bool ParseProblem(const std::string& filename, SOLVER& s)
   {
     mother_machine_cell_tracking_input i;
     const bool read_suc = read_input(filename, i);
     construct_tracking_problem(i, s);
     return read_suc;
   }

   template<typename SOLVER>
   bool ParseProblemMotherMachine(const std::string& filename, SOLVER& s)
   {
     mother_machine_cell_tracking_input i;
     const bool read_suc = read_input(filename, i);
     construct_tracking_problem(i, s);

      auto& cell_tracking_constructor = s.template GetProblemConstructor<0>();
      // add exit constraints based on positions of cell detection hypotheses
      // coordinates are in format (upper,lower) with upper < lower (hence coordinates are inverted!)
      //std::cout << "Exit constraints disabled\n";
      assert(i.cell_detection_stat.size() == i.cell_position.size());
      for(INDEX t=0; t<i.cell_position.size(); ++t) {
        assert(i.cell_detection_stat[t].size() == i.cell_position[t].size());
        // possibly better: if there exists cell strictly between i and j, then no exit constraint needs to be added
        for(INDEX d1=0; d1<i.cell_position[t].size(); ++d1) {
          for(INDEX d2=d1+1; d2<i.cell_position[t].size(); ++d2) {
            auto upper_bound_d1 = i.cell_position[t][d1][1];
            auto lower_bound_d1 = i.cell_position[t][d1][0];
            auto upper_bound_d2 = i.cell_position[t][d2][1];
            auto lower_bound_d2 = i.cell_position[t][d2][0];
            assert(upper_bound_d1 < lower_bound_d1);
            assert(upper_bound_d2 < lower_bound_d2);
            //std::cout << "t=" << t << ", H1=" << d1 << "(" << lower_bound_d1 << "," << upper_bound_d1 << "), H2=" << d2 << "(" << lower_bound_d2 << "," << upper_bound_d2 << ")\n";
            if(upper_bound_d1 > lower_bound_d2) {
              cell_tracking_constructor.add_exit_constraint(s.GetLP(), t,d1,d2);
            }
            if(upper_bound_d2 > lower_bound_d1) {
              cell_tracking_constructor.add_exit_constraint(s.GetLP(), t,d2,d1);
            }
          }
        } 
      }

      cell_tracking_constructor.order_factors(s.GetLP());

      return read_suc;
   }
} // end namespace cell_tracking_parser_mother_machine
*/

} // namespace LPMP

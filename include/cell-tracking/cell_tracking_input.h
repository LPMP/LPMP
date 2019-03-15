#ifndef LPMP_CELL_TRACKING_INPUT_H
#define LPMP_CELL_TRACKING_INPUT_H

#include <vector>
#include <cassert>
#include <utility>
#include <string>

namespace LPMP {

    struct cell_tracking_input {

        struct cell_detection {
            double detection_cost = 0.0, appearance_cost = 0.0, disappearance_cost = 0.0;
            std::size_t timestep, cell_number;
            std::size_t no_incoming_transition_edges = 0, no_incoming_division_edges = 0, no_outgoing_transition_edges = 0, no_outgoing_division_edges = 0;
            bool operator<(const cell_detection& o) const
            { 
                if(timestep != o.timestep) { return timestep < o.timestep; }
                else { return cell_number < o.cell_number; }
            }
        };
        // cell_detections must be ordere by timestep/cell_number (lexicographically)
        std::vector<cell_detection> cell_detections;
        std::size_t no_cells() const { return cell_detections.size(); }

        std::vector<std::size_t> timestep_to_first_cell_index; // given timestep, what is the first cell
        std::size_t no_timesteps() const { return timestep_to_first_cell_index.size(); }
        std::size_t no_cells(const std::size_t timestep) const {
            assert(timestep < no_timesteps());
            if(timestep < no_timesteps()-1) {
                return timestep_to_first_cell_index[timestep+1] - timestep_to_first_cell_index[timestep];
            } else {
                return no_cells() - timestep_to_first_cell_index[timestep];
            }
        }
        std::size_t cell_index(const std::size_t timestep, const std::size_t cell_number) const 
        {
            assert(timestep < no_timesteps());
            assert(cell_number < no_cells(timestep));
            return timestep_to_first_cell_index[timestep] + cell_number; 
        }

        struct cell_transition {
            std::size_t outgoing_cell;
            std::size_t incoming_cell;
            double cost; 
        };
        std::vector<cell_transition> cell_transitions;

        struct cell_division {
            std::size_t outgoing_cell;
            std::size_t incoming_cell_1, incoming_cell_2;
            double cost;
        };
        std::vector<cell_division> cell_divisions;

        std::vector<std::size_t> conflict_element_bounds;
        std::size_t no_conflicts() const { return conflict_element_bounds.size(); }

        std::vector<std::size_t> conflict_cells; 
        // for iterating over conflict sets
        using conflict_iterator_type = std::vector<std::size_t>::const_iterator;
        std::pair<conflict_iterator_type, conflict_iterator_type> get_conflict(const std::size_t conflict_no) const
        {
            assert(conflict_no < conflict_element_bounds.size());
            if(conflict_no < no_conflicts() -1) {
                return {conflict_cells.begin() + conflict_element_bounds[conflict_no], conflict_cells.begin() + conflict_element_bounds[conflict_no+1]};
            } else {
                return {conflict_cells.begin() + conflict_element_bounds[conflict_no], conflict_cells.end()};
            }
        }

        //std::vector<std::vector<detection_factor_stat>> cell_detection_stat; // detection cost, appearance cost, disappearance cost,  # incoming edges, # outgoing edges
        //std::vector<std::tuple<INDEX,INDEX,INDEX,INDEX,double>> mappings; // timestep outgoing, outgoing cell, timestep incoming, incoming cell, cost
        //std::vector<std::tuple<INDEX,INDEX,INDEX,INDEX,INDEX,INDEX,double>> divisions; // timestep outgoing, outgoing cell, incoming timestep 1, incoming cell 1, incoming timestep 2, incoming cell 2, cost
        //std::vector<std::vector<std::array<INDEX,2>>> conflicts; // timestep, {hyp_1, ..., hyp_n}
    };

    namespace cell_tracking_parser_2d {

        cell_tracking_input parse_file(const std::string& filename);
        cell_tracking_input parse_string(const std::string& filename);

    } 

} // namespace LPMP

#endif // LPMP_CELL_TRACKING_INPUT_H 

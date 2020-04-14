#pragma once

#include <vector>
#include <cassert>
#include <utility>
#include <string>
#include <numeric>

namespace LPMP {

    struct cell_tracking_instance {

        struct cell_detection {
            double detection_cost = 0.0, appearance_cost = 0.0, disappearance_cost = 0.0;
            std::size_t timestep = std::numeric_limits<std::size_t>::max(), cell_number = std::numeric_limits<std::size_t>::max();
            //std::size_t nr_incoming_transition_edges = 0, nr_incoming_division_edges = 0, nr_outgoing_transition_edges = 0, nr_outgoing_division_edges = 0;
            bool operator<(const cell_detection& o) const
            { 
                if(timestep != o.timestep) { return timestep < o.timestep; }
                else { return cell_number < o.cell_number; }
            }

            bool is_initial() const
            {
                return timestep == std::numeric_limits<std::size_t>::max() && cell_number == std::numeric_limits<std::size_t>::max();
            }
        };

        void add_cell_detection(const std::size_t timestep, const std::size_t hypothesis_id, const double cost);
        // cell_detections must be ordered by timestep/cell_number (lexicographically)
        std::vector<cell_detection> cell_detections;
        std::size_t nr_cells() const { return cell_detections.size(); }

        std::vector<std::size_t> timestep_to_first_cell_index; // given timestep, what is the first cell
        std::size_t nr_timesteps() const { return timestep_to_first_cell_index.size(); }
        std::size_t nr_cells(const std::size_t timestep) const {
            assert(timestep < nr_timesteps());
            if(timestep < nr_timesteps()-1) {
                return timestep_to_first_cell_index[timestep+1] - timestep_to_first_cell_index[timestep];
            } else {
                return nr_cells() - timestep_to_first_cell_index[timestep];
            }
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
        std::size_t nr_conflicts() const { return conflict_element_bounds.size(); }

        std::vector<std::size_t> conflict_cells; 
        // for iterating over conflict sets
        using conflict_iterator_type = std::vector<std::size_t>::const_iterator;
        std::pair<conflict_iterator_type, conflict_iterator_type> get_conflict(const std::size_t conflict_nr) const
        {
            assert(conflict_nr < conflict_element_bounds.size());
            if(conflict_nr < nr_conflicts() -1) {
                return {conflict_cells.begin() + conflict_element_bounds[conflict_nr], conflict_cells.begin() + conflict_element_bounds[conflict_nr+1]};
            } else {
                return {conflict_cells.begin() + conflict_element_bounds[conflict_nr], conflict_cells.end()};
            }
        }

        template<typename STREAM>
        void write_to_lp(STREAM& s) const;
    };

    namespace cell_tracking_parser_2d {

        cell_tracking_instance parse_file(const std::string& filename);
        cell_tracking_instance parse_string(const std::string& filename);

    }

    inline void cell_tracking_instance::add_cell_detection(const std::size_t timestep, const std::size_t hypothesis_id, const double cost)
    {
        if (hypothesis_id >= cell_detections.size())
            cell_detections.resize(hypothesis_id + 1);
        assert(cell_detections[hypothesis_id].is_initial());
        cell_detections[hypothesis_id].timestep = timestep;
        cell_detections[hypothesis_id].cell_number = hypothesis_id;
        cell_detections[hypothesis_id].detection_cost = cost;
    }

    template <typename STREAM>
    void cell_tracking_instance::write_to_lp(STREAM &s) const
    {
        auto cell_detection_identifier = [&](const cell_detection &c) -> std::string {
            assert(!c.is_initial());
            return std::string("det_" + std::to_string(c.timestep) + "_" + std::to_string(c.cell_number));
        };

        auto cell_appearance_identifier = [&](const cell_detection &c) -> std::string {
            return std::string("app_" + std::to_string(c.timestep) + "_" + std::to_string(c.cell_number));
        };

        auto cell_disappearance_identifier = [&](const cell_detection &c) -> std::string {
            return std::string("disapp_" + std::to_string(c.timestep) + "_" + std::to_string(c.cell_number));
        };

        auto cell_transition_identifier = [&](const cell_transition &t) -> std::string {
            const auto c_out = this->cell_detections[t.outgoing_cell];
            const auto c_in = this->cell_detections[t.incoming_cell];
            return std::string("trans_") + std::to_string(c_out.timestep) + "_" + std::to_string(c_out.cell_number) + "_to_" + std::to_string(c_in.timestep) + "_" + std::to_string(c_in.cell_number);
        };

        auto cell_division_identifier = [&](const cell_division &d) -> std::string {
            const auto c_out = this->cell_detections[d.outgoing_cell];
            const auto c_in_1 = this->cell_detections[d.incoming_cell_1];
            const auto c_in_2 = this->cell_detections[d.incoming_cell_2];
            return std::string("div_") + std::to_string(c_out.timestep) + "_" + std::to_string(c_out.cell_number) + "_to_" + std::to_string(c_in_1.timestep) + "_" + std::to_string(c_in_1.cell_number) + "_" + std::to_string(c_in_2.timestep) + "_" + std::to_string(c_in_2.cell_number);
        };

        s << "Minimize\n";
        for (const auto d : cell_detections)
        {
            if (!d.is_initial())
            {
                s << (d.detection_cost >= 0.0 ? "+ " : "") << d.detection_cost << " " << cell_detection_identifier(d) << "\n";
                s << (d.appearance_cost >= 0.0 ? "+ " : "") << d.appearance_cost << " " << cell_appearance_identifier(d) << "\n";
                s << (d.disappearance_cost >= 0.0 ? "+ " : "") << d.disappearance_cost << " " << cell_disappearance_identifier(d) << "\n";
            }
        }
        for (const auto &t : cell_transitions)
            s << (t.cost >= 0.0 ? "+ " : "") << t.cost << " " << cell_transition_identifier(t) << "\n";
        for (const auto &d : cell_divisions)
            s << (d.cost >= 0.0 ? "+ " : "") << d.cost << " " << cell_division_identifier(d) << "\n";

        s << "Subject To\n";
        std::vector<std::vector<std::size_t>> outgoing_cell_transitions(cell_detections.size());
        std::vector<std::vector<std::size_t>> incoming_cell_transitions(cell_detections.size());
        std::vector<std::vector<std::size_t>> outgoing_cell_divisions(cell_detections.size());
        std::vector<std::vector<std::size_t>> incoming_cell_divisions(cell_detections.size());
        for (std::size_t i = 0; i < cell_transitions.size(); ++i)
        {
            outgoing_cell_transitions[cell_transitions[i].outgoing_cell].push_back(i);
            outgoing_cell_transitions[cell_transitions[i].incoming_cell].push_back(i);
        }
        for (std::size_t i = 0; i < cell_divisions.size(); ++i)
        {
            outgoing_cell_divisions[cell_divisions[i].outgoing_cell].push_back(i);
            incoming_cell_divisions[cell_divisions[i].incoming_cell_1].push_back(i);
            incoming_cell_divisions[cell_divisions[i].incoming_cell_2].push_back(i);
        }
        // incoming flow conservation
        for (std::size_t i = 0; i < cell_detections.size(); ++i)
        {
            if(!cell_detections[i].is_initial())
            {
                s << " - " << cell_detection_identifier(cell_detections[i]) << " + " << cell_appearance_identifier(cell_detections[i]);
                for (const std::size_t t : incoming_cell_transitions[i])
                    s << " + " << cell_transition_identifier(cell_transitions[t]);
                for (const std::size_t d : incoming_cell_divisions[i])
                    s << " + " << cell_division_identifier(cell_divisions[d]);
                s << " = 0\n";
            }
        }

        // outgoing flow conservation
        for (std::size_t i = 0; i < cell_detections.size(); ++i)
        {
            if(!cell_detections[i].is_initial())
            {
                s << " - " << cell_detection_identifier(cell_detections[i]) << " + " << cell_disappearance_identifier(cell_detections[i]);
                for (const std::size_t t : outgoing_cell_transitions[i])
                    s << " + " << cell_transition_identifier(cell_transitions[t]);
                for (const std::size_t d : outgoing_cell_divisions[i])
                    s << " + " << cell_division_identifier(cell_divisions[d]);
                s << " = 0\n";
            }
        }

        // exclusion constraints
        for (std::size_t conflict_nr = 0; conflict_nr < nr_conflicts(); ++conflict_nr)
        {
            auto [conflict_begin, conflict_end] = get_conflict(conflict_nr);
            for (auto it = conflict_begin; it != conflict_end; ++it)
                s << " + " << cell_detection_identifier(cell_detections[*it]);
            s << " <= 1\n";
        }

        s << "Bounds\nBinaries\n";
        for (const auto d : cell_detections)
        {
            if(!d.is_initial())
            {
                s << cell_detection_identifier(d) << "\n";
                s << cell_appearance_identifier(d) << "\n";
                s << cell_disappearance_identifier(d) << "\n";
            }
        }
        for (const auto &t : cell_transitions)
            s << cell_transition_identifier(t) << "\n";
        for (const auto &d : cell_divisions)
            s << cell_division_identifier(d) << "\n";

        s << "End\n";
        }

} // namespace LPMP
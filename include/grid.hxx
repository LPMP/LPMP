#ifndef LPMP_GRID_HXX
#define LPMP_GRID_HXX

#include <array>
#include <vector>
#include <cassert>
#include <type_traits>

namespace LPMP {

    struct grid {

        std::array<std::size_t,2> node_to_grid(const std::size_t i) const 
        {
            assert(i < size_x*size_y);
            const std::size_t x = i%size_x;
            const std::size_t y = i/size_x;
            return {x,y}; 
        }

        // chain number, position in chain
        const std::array<std::size_t,2> edge_to_chain(const std::size_t node_1, const std::size_t node_2) const
        {
            if(size_x == 1 || size_y == 1) {
                assert(node_2 == node_1 + 1 && node_2 < size_x*size_y);
                return {0,node_1}; 
            } else if(node_2 == node_1 + 1) { // Horizontally arranged
                const auto coord_1 = node_to_grid(node_1);
                const auto coord_2 = node_to_grid(node_2);
                assert(coord_2[0] == coord_1[0] + 1);
                assert(coord_2[1] == coord_1[1]);
                return {coord_1[1], coord_1[0]};
            } else if(node_2 == node_1 + size_x) { // Vertically arranged
                const auto coord_1 = node_to_grid(node_1);
                const auto coord_2 = node_to_grid(node_2);
                assert(coord_2[1] == coord_1[1] + 1);
                assert(coord_1[0] == coord_2[0]);
                return {coord_1[0] + size_y, coord_1[1]}; 
            } else {
                throw std::runtime_error("nodes must be neighboring in grid graph");
            }
        }
        
        const std::size_t number_of_chains() const 
        { 
            if(size_x == 1 || size_y == 1) { return 1; }
            else { return size_x + size_y; } 
        }

        std::vector<std::size_t> chain(const std::size_t chain_no) const
        {
            assert(chain_no < number_of_chains());
            if(size_x == 1 || size_y == 1) {
                std::vector<std::size_t> c(size_x*size_y);
                std::iota(c.begin(), c.end(), 0);
                return c;
            } else {
                if(chain_no < size_y) { return horizontal_chain(chain_no); }
                else { return vertical_chain(chain_no - size_y); }
            }
        }

        std::vector<std::size_t> horizontal_chain(const std::size_t j) const 
        {
            assert(j < size_y);
            std::vector<std::size_t> nodes;
            nodes.reserve(size_x);
            for(std::size_t i=0; i<size_x; i++) {
                nodes.push_back(j*size_x + i);
            }
            return nodes;
        }

        std::vector<std::size_t> vertical_chain(const std::size_t i) const 
        {
            assert(i < size_x);
            std::vector<std::size_t> nodes;
            nodes.reserve(size_y);
            for(std::size_t j=0; j<size_y; j++) {
                nodes.push_back(j*size_x + i);
            }
            return nodes;
        }

        const std::size_t size_x;
        const std::size_t size_y;
    };

    grid recognize_grid(const std::vector<std::array<std::size_t,2>>& edges)
    {
        const std::size_t x_edge_distance = 1;
        std::size_t y_edge_distance = 0;
        std::size_t last_horizontal_node = 0;
        std::size_t horizontal_length = 0;
        std::size_t last_vertical_node = 0;
        std::size_t vertical_length = 0;

        for(const auto& e : edges) {
            assert(e[1] > e[0]);
            const std::size_t current_edge_distance = e[1] - e[0];
            if (y_edge_distance == 0 && current_edge_distance != x_edge_distance) {
                y_edge_distance = current_edge_distance;
            } 
            if (current_edge_distance != x_edge_distance && y_edge_distance != current_edge_distance) {
                throw std::runtime_error("could not recognize graph as grid");
            }

            if (e[0] == last_horizontal_node && current_edge_distance == x_edge_distance) {
                horizontal_length++;
                last_horizontal_node = e[1];
            }

            if (e[0] == last_vertical_node && current_edge_distance == y_edge_distance) {
                vertical_length++;
                last_vertical_node = e[1];
            }
        }
        const std::size_t grid_size_x = horizontal_length + 1;
        const std::size_t grid_size_y = vertical_length + 1;

        return grid({grid_size_x, grid_size_y});
    }

} // namespace LPMP

#endif // LPMP_GRID_HXX


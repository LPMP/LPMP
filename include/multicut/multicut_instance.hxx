#ifndef LPMP_MULTICUT_INSTANCE
#define LPMP_MULTICUT_INSTANCE

namespace LPMP {

   struct multicut_instance {
      struct weighted_edge : public std::array<std::size_t,2> { 
         weighted_edge()
            : std::array<std::size_t,2>({0,0}), cost(0.0) {}

         weighted_edge(const std::size_t i, const std::size_t j, const double c)
            : std::array<std::size_t,2>({i,j}), cost(c) {}

         double cost;
      };

      void add_edge(const std::size_t i, const std::size_t j, const double cost)
      {
         assert(i != j);
         no_nodes = std::max({no_nodes,i+1,j+1});
         edges.push_back(weighted_edge(i,j,cost));
      }

      std::size_t no_nodes = 0;
      std::vector<weighted_edge> edges; 
   };

} // namespace LPMP

#endif // LPMP_MULTICUT_INSTANCE

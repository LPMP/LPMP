#ifndef LPMP_BINARY_MRF_INSTANCE_HXX
#define LPMP_BINARY_MRF_INSTANCE_HXX

#include <vector>
#include <array>
#include <tuple>
#include <cassert>

namespace LPMP {

struct binary_MRF_instance {
    struct binary_pairwise_potential {
        binary_pairwise_potential(const std::size_t _i, const std::size_t _j, std::array<std::array<double,2>,2> _cost = {{{0.0,0.0},{0.0,0.0}}})
            : i(_i), j(_j), cost(_cost) {}
        binary_pairwise_potential() = default;

        using msg_type = std::array<double,2>;

        std::tuple<double, msg_type, msg_type> make_potts() const
        {
            const auto c = 0.5*(cost[0][1] + cost[1][0] - cost[0][0] - cost[1][1]);
            msg_type msg_1 {{0.0, cost[1][0] - c - cost[0][0]}};
            msg_type msg_2 {{cost[0][0], cost[0][1] - c}};

            return {c, msg_1, msg_2};
        }

        void reparametrize(const msg_type& msg_1, const msg_type& msg_2)
        {
            cost[0][0] -= msg_1[0] + msg_2[0];
            cost[0][1] -= msg_1[0] + msg_2[1];
            cost[1][0] -= msg_1[1] + msg_2[0];
            cost[1][1] -= msg_1[1] + msg_2[1];
        }

        bool has_same_support(const binary_pairwise_potential& o) const
        {
            return i == o.i && j == o.j;
        }

        void add_cost(const binary_pairwise_potential& o)
        {
            cost[0][0] += o.cost[0][0];
            cost[0][1] += o.cost[0][1];
            cost[1][0] += o.cost[1][0];
            cost[1][1] += o.cost[1][1];
        }

        void transpose()
        {
           std::swap(i,j);
           std::swap(cost[0][1], cost[1][0]);
        }

        bool operator<(const binary_pairwise_potential& o) const
        {
           std::array<std::size_t,2> idx1 = {i,j};
           std::array<std::size_t,2> idx2 = {o.i,o.j};
           return std::lexicographical_compare(idx1.begin(), idx1.end(), idx2.begin(), idx2.end());
        }

        std::size_t i = 0;
        std::size_t j = 0;
        std::array<std::array<double,2>,2> cost = {{ {0.0, 0.0}, {0.0, 0.0} }};
    };

    double constant = 0.0;
    std::vector<std::array<double,2>> unaries;
    std::vector<binary_pairwise_potential> pairwise_potentials;

    using labeling = std::vector<unsigned char>;

    double evaluate(const labeling& l) const
    {
        assert(l.size() == unaries.size());
        double cost = constant;
        for(std::size_t i=0; i<unaries.size(); ++i) {
            assert(l[i] == 0 || l[i] == 1);
            cost += unaries[i][l[i]];
        }

        for(const auto pot : pairwise_potentials) {
            cost += pot.cost[ l[pot.i] ][ l[pot.j] ];
        }
        
        return cost;
    }

    // write file in uai format
    template<typename STREAM>
    void write_uai(STREAM& s) const
    {
        // preamble
        s << "MARKOV\n";
        s << unaries.size() << "\n";
        // cadinalities
        for(std::size_t i=0; i<unaries.size(); ++i) {
            s << 2;
            if(i+1<unaries.size()) s << " ";
            else s << "\n";
        }
        // number of potentials
        s << unaries.size() + pairwise_potentials.size() << "\n";
        // variables of potentials
        for(std::size_t i=0; i<unaries.size(); ++i) {
            s << 1 << " " << i << "\n";
        }
        for(const auto& p : pairwise_potentials) {
            s << 2 << " " << p.i << " " << p.j << "\n";
        }

        // function tables
        for(const auto& u : unaries) {
            s << 2 << " " << u[0] << " " << u[1] << "\n";
        }
        for(const auto& p : pairwise_potentials) {
            s << 4 << " " << p.cost[0][0] << " " << p.cost[0][1] << " " << p.cost[1][0] << " " << p.cost[1][1] << "\n";
        } 
    }
};

struct binary_Potts_instance {
    struct weighted_edge : public std::array<std::size_t,2> { 
        weighted_edge(const std::size_t i, const std::size_t j, const double _cost) : std::array<std::size_t,2>({i,j}), cost(_cost) {}
        double cost; 
    };

    double constant = 0.0;
    std::vector<std::array<double,2>> unaries;
    std::vector<weighted_edge> pairwise_potentials;

    using labeling = std::vector<unsigned char>;

    double evaluate(const labeling& l) const
    {
        assert(l.size() == unaries.size());
        double cost = constant;
        for(std::size_t i=0; i<unaries.size(); ++i) {
            assert(l[i] == 0 || l[i] == 1);
            cost += unaries[i][l[i]];
        }

        for(const auto pot : pairwise_potentials) {
            if(l[pot[0]] != l[pot[1]]) {
                cost += pot.cost;
            }
        }
        
        return cost;
    }

    // write file in uai format
    template<typename STREAM>
    void write_uai(STREAM& s) const
    {
        // preamble
        s << "MARKOV\n";
        s << unaries.size() << "\n";
        for(std::size_t i=0; i<unaries.size(); ++i) {
            s << 2;
            if(i+1<unaries.size()) s << " ";
            else s << "\n";
        }
        s << unaries.size() + pairwise_potentials.size() << "\n";
        for(std::size_t i=0; i<unaries.size(); ++i) {
            s << 1 << " " << i << "\n";
        }
        for(const auto& p : pairwise_potentials) {
            s << 2 << " " << p[0] << " " << p[1] << "\n";
        }

        // function tables
        for(const auto& u : unaries) {
            s << 2 << " " << u[0] << " " << u[1] << "\n";
        }
        for(const auto& p : pairwise_potentials) {
            s << 4 << " " << 0.0 << " " << p.cost << " " << p.cost << " " << 0.0 << "\n";
        } 
    }
};

binary_Potts_instance transform_binary_MRF_to_Potts(const binary_MRF_instance& binary_MRF)
{
    binary_Potts_instance potts;
    potts.unaries = binary_MRF.unaries;

    for(const auto& pairwise_pot : binary_MRF.pairwise_potentials) {
        // see Schraudolph2009, eq. (3)
        const std::size_t i = pairwise_pot.i;
        const std::size_t j = pairwise_pot.j;
        const auto& p = pairwise_pot.cost;
        const auto c = 0.5*(p[0][1] + p[1][0] - p[0][0] - p[1][1]);
        const auto delta_i = p[1][0] - p[0][0] - c;
        const auto delta_j = p[0][1] - p[0][0] - c;

        potts.unaries[i][0] -= delta_i;
        potts.unaries[j][0] -= delta_j;
        binary_Potts_instance::weighted_edge e(i,j,c);
        potts.pairwise_potentials.push_back(e);

        //potts.constant += // TODO
    }

    return potts;
}

} // namespace LPMP

#endif // LPMP_BINARY_MRF_INSTANCE_HXX

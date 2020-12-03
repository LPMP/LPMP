#include "test.h"
#include "polynomial_functions.h"
#include <algorithm>
#include <random>

using namespace LPMP;

int main()
{
    // cubic function roots
    {
        cubic_function c(1,4,-8,9);
        auto roots = c.roots();
        test(std::abs(roots[0] - -5.685508364373402) < 1e-8);
        test(roots[1] == std::numeric_limits<double>::infinity());
        test(roots[2] == std::numeric_limits<double>::infinity()); 
    }

    {
        cubic_function c(1,-4,1,6);
        auto roots = c.roots();
        std::sort(roots.begin(), roots.end());
        test(std::abs(roots[0] - -1) < 1e-8);
        test(std::abs(roots[1] - 2) < 1e-8);
        test(std::abs(roots[2] - 3) < 1e-8);
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dis(0.0, 1.0);
    for(std::size_t i=0; i<100; ++i) {
        cubic_function f(dis(gen), dis(gen), dis(gen), dis(gen));
        const auto roots = f.roots();
        for(const double x : roots) {
            if(x != std::numeric_limits<double>::infinity()) {
                test(std::abs(f.evaluate(x)) <= 1e-8);
            }
        } 
    }
}

#include <array>
#include <cmath>

const double PI = 3.141592653589793238463;

namespace LPMP {

    class cubic_function {
        public:
            cubic_function(const double cubic_term, const double quadratic_term, const double linear_term, const double constant)
                : A(cubic_term),
                B(quadratic_term),
                C(linear_term),
                D(constant)
        {}

            double evaluate(const double x) const
            {
                return A*std::pow(x,3) + B*std::pow(x,2) + C*x + D;
            }

            std::array<double,3> roots() const
            {
                // bring to standard formula
                const double a = B/A;
                const double b = C/A;
                const double c = D/A;

                // eliminate quadratic term by substituting x <- x - a/3,
                // obtain x^3 + px + q
                const double p = b-std::pow(a,2)/3.0;
                const double q = 2*std::pow(a,3)/27.0 - a*b/3.0 + c;

                const double discriminant = std::pow(q/2.0,2) + std::pow(p/3.0,3);

                if(discriminant > 0.0) {
                    const double u = std::cbrt(-q/2.0 + std::sqrt(discriminant));
                    const double v = std::cbrt(-q/2.0 - std::sqrt(discriminant));
                    return {u+v - B/(3*A), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};;
                } else if(discriminant < 0.0) {
                    const double root_1 = std::sqrt(-4.0/3.0*p) * (1.0/3.0 * std::acos(-q/2.0 * std::sqrt(-27/std::pow(p,3))) + PI/3.0 ) - B/(3.0*A);
                    const double root_2 = std::sqrt(-4.0/3.0*p) * (1.0/3.0 * std::acos(-q/2.0 * std::sqrt(-27/std::pow(p,3))) ) - B/(3.0*A);
                    const double root_3 = std::sqrt(-4.0/3.0*p) * (1.0/3.0 * std::acos(-q/2.0 * std::sqrt(-27/std::pow(p,3))) - PI/3.0 ) - B/(3.0*A);
                    return {root_1, root_2, root_3};
                } else {
                    assert(discriminant == 0.0);
                    if(p == 0) {
                        return {-B/(3.0*A), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
                    } else {
                        return {3*q/p - B/(3*A), -3.0*q/(2.0*p) - B/(3.0*A), std::numeric_limits<double>::infinity()};
                    }
                }
            }

        private:
                // ax^3 + bx^2 + cx + d
                const double A,B,C,D; 
    };

    struct quartic_function {
        public:
            quartic_function(const double quartic_term, const double cubic_term, const double quadratic_term, const double linear_term, const double constant)
            : A(quartic_term),
            B(cubic_term),
            C(quadratic_term),
            D(linear_term),
            E(constant)
        {}

        double evaluate(const double x) const
        {
            return A*std::pow(x,4) + B*std::pow(x,3) + C*std::pow(x,2) + D*x + E;
        }

        double evaluate_derivative(const double x) const
        {
            cubic_function c(4*A, 3*B, 2*C, D);
            return c.evaluate(x); 
        }

        std::array<double,2> minima() const
        {
            cubic_function c(4*A, 3*B, 2*C, D);
            auto roots = c.roots();

            std::array<double,2> minima;
            std::size_t minima_idx = 0;
            for(std::size_t i=0; i<roots.size(); ++i) {
                const double x = roots[i];
                if(c.evaluate(x) >= 0.0) {
                    minima[minima_idx++] = x; 
                }
            }

            assert(minima_idx > 0);
            return minima;
        }

        private:
            const double A,B,C,D,E;
    };
}

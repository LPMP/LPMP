#include "test.h"
#include "vector.hxx"
#include <random>

using namespace LPMP;

template<typename VECTOR>
void test_vector_minima(const std::size_t n, std::random_device& rd)
{
    std::mt19937 gen{rd()};
    std::normal_distribution<> dist{5,2};

    VECTOR v(n);
    for(std::size_t i=0; i<n; ++i) {
        v[i] = dist(gen);
    }

    test(v.min() == *std::min_element(v.begin(), v.end())); 

    test(v.min_except(0) == *std::min_element(v.begin()+1, v.end()));
    for(std::size_t i=1; i<n-1; ++i) {
        const auto min_except = std::min( *std::min_element(v.begin(), v.begin() + i), *std::min_element(v.begin()+i+1, v.end()));
        test(v.min_except(i) == min_except);
    }
    test(v.min_except(n-1) == *std::min_element(v.begin(), v.begin()+n-1));

    const auto two_min = v.two_min();
    std::sort(v.begin(), v.end());
    test(v[0] == two_min[0]);
    test(v[1] == two_min[1]);
}

int main() {

  { // vector minimum

      LPMP::vector<REAL> v(5);
    v[0] = -1.0;
    v[1] = 0.0;
    v[2] = 1.0;
    v[3] = 2.0;
    v[4] = 3.0;

    test(v.min() == -1.0); 
    std::cout << v;

    test(v.min_except(0) == 0.0);
    test(v.min_except(1) == -1.0);

    auto two_min = v.two_min();
    test(two_min[0] == -1.0);
    test(two_min[1] == 0.0);
  }

  { // random vector tests
    std::random_device rd{};

    for(std::size_t n=2; n<100; ++n) {
        test_vector_minima<LPMP::vector<double>>(n, rd);
        test_vector_minima<min_vector<double>>(n, rd);
    }
  }

  { // matrix minima
    matrix<REAL> m(5,6);
    m(0,0) = -2.0; m(0,1) = +0.0; m(0,2) = +2.0; m(0,3) = -0.5; m(0,4) = +0.0; m(0,5) = +0.5;
    m(1,0) = -1.0; m(1,1) = +0.0; m(1,2) = +1.0; m(1,3) = -0.5; m(1,4) = +0.0; m(1,5) = +0.5;
    m(2,0) = -0.0; m(2,1) = -4.0; m(2,2) = +0.5; m(2,3) = -0.5; m(2,4) = +0.0; m(2,5) = +0.5;
    m(3,0) = +1.0; m(3,1) = +0.0; m(3,2) = -1.0; m(3,3) = -0.5; m(3,4) = +0.0; m(3,5) = +0.5;
    m(4,0) = +2.0; m(4,1) = +0.0; m(4,2) = -2.0; m(4,3) = -0.5; m(4,4) = +0.0; m(4,5) = +0.5;

    std::cout << m;

    { // matrix columnwise minimum
      auto min_col = m.min2(); 
      test(min_col.size() == 6);

      test(min_col[0] == -2.0);
      test(min_col[1] == -4.0);
      test(min_col[2] == -2.0); 
      test(min_col[3] == -0.5); 
      test(min_col[4] == +0.0); 
      test(min_col[5] == +0.5); 
    }

    { // matrix rowwise minimum
      auto min_row = m.min1();
      test(min_row.size() == 5);

      test(min_row[0] == -2.0);
      test(min_row[1] == -1.0);
      test(min_row[2] == -4.0);
      test(min_row[3] == -1.0);
      test(min_row[4] == -2.0); 
    }
  } 
}


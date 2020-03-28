#ifndef LPMP_TEST_H
#define LPMP_TEST_H

#include <stdexcept>
#include <string>
#include <iostream>

const std::string red("\033[0;31m");
const std::string green("\033[1;32m");
const std::string yellow("\033[1;33m");
const std::string cyan("\033[0;36m");
const std::string magenta("\033[0;35m");
const std::string reset("\033[0m");

void test(const double given_value, const double expected_value, double percentage_tolerance = 10)
{
  const double eps = 1e-8;

  if (percentage_tolerance < eps)
    percentage_tolerance = eps; // To handle numerical errors

  auto percentage_difference = 100 * std::abs(given_value - expected_value) / given_value;
  if (given_value == expected_value)
    percentage_difference = 0; // To handle infinities

  if (percentage_difference <= percentage_tolerance)
    std::cout<<green<<"Test Passed With Tolerance of: "<<cyan<<percentage_difference<<"% "<<green<<", Within expected tolerance: "<<magenta<<percentage_tolerance<<"% "<<reset<<std::endl;
  else
  {
    std::cout<<red<<"Test Failed With Tolerance of: "<<cyan<<percentage_difference<<"% "<<red<<", Outside expected tolerance: "<<magenta<<percentage_tolerance<<"% "<<reset<<std::endl;
    throw std::runtime_error("Test failed.");
  }
}

#endif // LPMP_TEST_H

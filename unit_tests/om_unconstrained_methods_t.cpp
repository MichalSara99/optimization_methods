#define BOOST_TEST_MODULE

#ifdef _DEBUG_
#include <iostream>
#endif

#include "unconstrained_methods/om_unconstrained_methods.hpp"
#include "utilities/om_types.hpp"
#include "utilities/om_utilities.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

using om_types::f_scalar_t;
using om_unconstrained_methods::minimize;
using om_utilities::range;

BOOST_AUTO_TEST_SUITE(minimize_one_dimensional)

BOOST_AUTO_TEST_CASE(brent_default) {

  auto minF = [](double x) -> double {
    return (x * x + 2.0 * std::exp(-1.0 * x));
  };
  auto r = range<>{0.0, 2.0};
  std::size_t max_iters = 1000;
  auto tol = 10e-5;
  auto g_min = minimize<>((f_scalar_t<double>)minF, r, tol, max_iters);

  // check if we didnt reech the number of iteration:
  BOOST_CHECK_LT(std::get<3>(g_min), max_iters);
  // minF is obviously always > 0 so lets check the minimum value
  BOOST_CHECK_GT(std::get<1>(g_min), tol);

#define _DEBUG_

#ifdef _DEBUG_
  std::cout << "================== Brent Method (Default) ==============\n";
  std::cout << "minimiser: " << std::get<0>(g_min) << "\n"
            << "fun value: " << std::get<1>(g_min) << "\n"
            << "fun evaluations: " << std::get<2>(g_min) << "\n"
            << "iterations: " << std::get<3>(g_min) << "\n";
#endif
}

BOOST_AUTO_TEST_CASE(brent_float) {

  auto minF = [](float x) { return (x * x + 2.0 * std::exp(-1.0 * x)); };
  auto r = range<float>{0.0, 2.0};
  std::size_t max_iters = 1000;
  auto tol = 10e-5;
  auto g_min = minimize<float>((f_scalar_t<float>)minF, r, tol, max_iters);

  // check if we didnt reech the number of iteration:
  BOOST_CHECK_LT(std::get<3>(g_min), max_iters);
  // minF is obviously always > 0 so lets check the minimum value
  BOOST_CHECK_GT(std::get<1>(g_min), tol);

#define _DEBUG_

#ifdef _DEBUG_
  std::cout << "================== Brent Method (Default) ==============\n";
  std::cout << "minimiser: " << std::get<0>(g_min) << "\n"
            << "fun value: " << std::get<1>(g_min) << "\n"
            << "fun evaluations: " << std::get<2>(g_min) << "\n"
            << "iterations: " << std::get<3>(g_min) << "\n";
#endif
}

BOOST_AUTO_TEST_CASE(brent_default_1) {

  auto minF = [](double x) { return (x * std::cos(x)); };
  const double pi{3.14159265359};
  auto r = range<>{0.0, pi / 2.0};
  auto tol = 10e-5;
  std::size_t max_iters = 1000;
  auto g_min = minimize<>((f_scalar_t<double>)minF, r, tol, max_iters);

  // check if we didnt reech the number of iteration:
  BOOST_CHECK_LT(std::get<3>(g_min), max_iters);
  // the minimum is at x = 0 and minF = 0 there
  // lets check those
  BOOST_CHECK_LT(std::abs(std::get<0>(g_min)), tol);
  BOOST_CHECK_LT(std::abs(std::get<1>(g_min)), tol);

#ifdef _DEBUG_
  std::cout << "================== Brent Method (Default) ==============\n";
  std::cout << "minimiser: " << std::get<0>(g_min) << "\n"
            << "fun value: " << std::get<1>(g_min) << "\n"
            << "fun evaluations: " << std::get<2>(g_min) << "\n"
            << "iterations: " << std::get<3>(g_min) << "\n";

#endif
}

BOOST_AUTO_TEST_SUITE_END()

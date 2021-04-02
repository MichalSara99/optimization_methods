#define BOOST_TEST_MODULE

#ifdef _DEBUG_
#include <iostream>
#endif

#include "unconstrained_methods/multi_dim/steepest_descent/om_steepest_descent.hpp"
#include "unconstrained_methods/multi_dim/test_functions/om_raos_collection.hpp"
#include "unconstrained_methods/multi_dim/test_functions/om_test_helpers.hpp"
#include "unconstrained_methods/one_dim/om_line_methods.hpp"
#include "utilities/om_types.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

using namespace om_unconstrained_methods::om_line_methods;
using namespace om_unconstrained_methods::om_steepest_descent;
using namespace om_test_helpers;

BOOST_AUTO_TEST_SUITE(steepest_descent_test)

BOOST_AUTO_TEST_CASE(steepest_descent_with_gs_default) {

  auto helper = create_rao_test_collection<double>();

#define _DEBUG_

#ifdef _DEBUG_
  std::cout << "\n\n";
  std::cout << "Running Rao tests.... \n";
#endif
  auto golden_section =
      golden_section_method<double>{range<double>{-1.0, 1.0}, 1.0e-7};
  std::size_t multi_max_iters = 10000;
  auto multi_dim_method = steepest_descent_method<double>{
      golden_section, multi_max_iters, 1.0e-7, 1.0e-7, 1.0e-7};
  double tol = 1.0e-3;
  for (auto const &h : helper) {
    auto result = multi_dim_method.minimize(h->objective, h->guess);
    const auto &expected_min = h->minimizer;
    const auto &found_min = std::get<0>(result);
    BOOST_CHECK_LE((expected_min - found_min).norm(), tol);
    BOOST_CHECK_LT(std::get<2>(result), multi_max_iters);

#ifdef _DEBUG_
    std::cout << "FUNCTION NAME: " << h->name << "\n";
    std::cout << "EXPECTED MINIMISER: \n" << h->minimizer << "\n";
    std::cout << "FOUND MINIMISER: \n" << std::get<0>(result) << "\n";
    std::cout << "FUN VALUE: " << std::get<1>(result) << "\n";
    std::cout << "ITERATIONS: " << std::get<2>(result) << "\n";
    std::cout << "---------------------------------------------------\n";
#endif
  }
}

BOOST_AUTO_TEST_CASE(steepest_descent_with_f_default) {

  auto helper = create_rao_test_collection<double>();

#define _DEBUG_

#ifdef _DEBUG_
  std::cout << "\n\n";
  std::cout << "Running Rao tests.... \n";
#endif
  auto fibonacci = fibonacci_method<double>{range<double>{-1.0, 1.0}, 1.0e-7};
  std::size_t multi_max_iters = 10000;
  auto multi_dim_method = steepest_descent_method<double>{
      fibonacci, multi_max_iters, 1.0e-7, 1.0e-7, 1.0e-7};
  double tol = 1.0e-3;
  for (auto const &h : helper) {
    auto result = multi_dim_method.minimize(h->objective, h->guess);
    const auto &expected_min = h->minimizer;
    const auto &found_min = std::get<0>(result);
    BOOST_CHECK_LE((expected_min - found_min).norm(), tol);
    BOOST_CHECK_LT(std::get<2>(result), multi_max_iters);

#ifdef _DEBUG_
    std::cout << "FUNCTION NAME: " << h->name << "\n";
    std::cout << "EXPECTED MINIMISER: \n" << h->minimizer << "\n";
    std::cout << "FOUND MINIMISER: \n" << std::get<0>(result) << "\n";
    std::cout << "FUN VALUE: " << std::get<1>(result) << "\n";
    std::cout << "ITERATIONS: " << std::get<2>(result) << "\n";
    std::cout << "---------------------------------------------------\n";
#endif
  }
}

BOOST_AUTO_TEST_CASE(steepest_descent_with_b_default) {

  auto helper = create_rao_test_collection<double>();

#define _DEBUG_

#ifdef _DEBUG_
  std::cout << "\n\n";
  std::cout << "Running Rao tests.... \n";
#endif
  auto brent = brent_method<double>{range<double>{-1.0, 1.0}, 1.0e-7};
  std::size_t multi_max_iters = 10000;
  auto multi_dim_method = steepest_descent_method<double>{
      brent, multi_max_iters, 1.0e-7, 1.0e-7, 1.0e-7};
  double tol = 1.0e-1;
  for (auto const &h : helper) {
    auto result = multi_dim_method.minimize(h->objective, h->guess);
    const auto &expected_min = h->minimizer;
    const auto &found_min = std::get<0>(result);
    BOOST_CHECK_LE((expected_min - found_min).norm(), tol);
    BOOST_CHECK_LT(std::get<2>(result), multi_max_iters);

#ifdef _DEBUG_
    std::cout << "FUNCTION NAME: " << h->name << "\n";
    std::cout << "EXPECTED MINIMISER: \n" << h->minimizer << "\n";
    std::cout << "FOUND MINIMISER: \n" << std::get<0>(result) << "\n";
    std::cout << "FUN VALUE: " << std::get<1>(result) << "\n";
    std::cout << "ITERATIONS: " << std::get<2>(result) << "\n";
    std::cout << "---------------------------------------------------\n";
#endif
  }
}

BOOST_AUTO_TEST_SUITE_END()

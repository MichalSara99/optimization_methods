#define BOOST_TEST_MODULE

#ifdef _DEBUG_
#include <iostream>
#endif

#include "unconstrained_methods/multi_dim/test_functions/om_raos_collection.hpp"
#include "unconstrained_methods/multi_dim/test_functions/om_test_helpers.hpp"
#include "unconstrained_methods/multi_dim/zero_order/om_nelder_mead.hpp"
#include "unconstrained_methods/one_dim/om_line_methods.hpp"
#include "utilities/om_types.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

using namespace om_unconstrained_methods::om_line_methods;
using namespace om_unconstrained_methods::om_zero_order;
using namespace om_test_helpers;

BOOST_AUTO_TEST_SUITE(nelder_mead)

BOOST_AUTO_TEST_CASE(nelder_mead_default) {

  auto helper = create_rao_test_collection<double>();

#define _DEBUG_

#ifdef _DEBUG_
  std::cout << "\n\n";
  std::cout << "Running Rao tests.... \n";
#endif

  std::size_t multi_max_iters = 50000;
  auto multi_dim_method = nelder_mead_method<double>{multi_max_iters, 10e-8};
  double tol = 1.0e-2;
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

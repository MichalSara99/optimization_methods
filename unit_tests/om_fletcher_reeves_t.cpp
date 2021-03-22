#define BOOST_TEST_MODULE

#ifdef _DEBUG_
#include <iostream>
#endif

#include "unconstrained_methods/multi_dim/conjugate_gradient/om_fletcher_reeves.hpp"
#include "unconstrained_methods/multi_dim/test_functions/om_raos_collection.hpp"
#include "unconstrained_methods/multi_dim/test_functions/om_test_helpers.hpp"
#include "unconstrained_methods/one_dim/om_line_methods.hpp"
#include "utilities/om_types.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

using namespace om_unconstrained_methods::om_line_methods;
using namespace om_unconstrained_methods::om_conjugate_gradient;
using namespace om_test_helpers;

BOOST_AUTO_TEST_SUITE()

BOOST_AUTO_TEST_CASE(fletcher_reeves_with_gs_default) {

  auto helper = create_rao_test_collection<double>();

#ifdef _DEBUG_
  std::cout << "\n\n";
  std::cout << "Running Rao tests.... \n";
#endif
  auto golden_section = golden_section_method<double>{range<double>{0.0, 1.0}};
  auto multi_dim_method = fletcher_reeves_method<double>{
      golden_section, 300, 1.0e-4, 1.0e-4, 1.0e-4};
  double tol = 1.0e-4;
  for (auto const &h : helper) {
    auto result = multi_dim_method.minimize(h->objective, h->guess);

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

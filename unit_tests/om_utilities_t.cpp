#define BOOST_TEST_MODULE

#include "utilities/om_types.hpp"
#include "utilities/om_utilities.hpp"
#include <iostream>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

using namespace om_utilities;
using om_types::vector_arg_t;

BOOST_AUTO_TEST_SUITE()

BOOST_AUTO_TEST_CASE(rand_vector_normal) {

  auto initGuess = Eigen::Vector3d(1.0, 2.0, 3.5);
  auto container = random_vectors_from_guess<>()(4, initGuess);
  BOOST_CHECK(!container.empty());
  // max coefficient times 3 times standard error:
  double tol = 3.0 * initGuess.maxCoeff();

  for (auto const &v : container) {
    for (std::size_t r = 0; r < v.rows(); ++r) {
      BOOST_CHECK_LE(std::abs(v(r)), tol);
    }
  }
}

BOOST_AUTO_TEST_CASE(rand_vectors_specific_dist) {

  auto initGuess = Eigen::Vector3d(1.0, 2.0, 3.5);
  auto container =
      random_vectors_from_guess<double, std::uniform_real_distribution>()(
          4, initGuess);
  auto tol = initGuess.maxCoeff() * (1.0 + 1. / sqrt(12.0));

  for (auto const &v : container) {
    for (std::size_t r = 0; r < v.rows(); ++r) {
      BOOST_CHECK_LE(v(r), tol);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

void rand_vectors_specific_dist_sort() {
  std::cout << "Initial guess: \n";
  auto initGuess = Eigen::Vector3d(1.0, 2.0, 3.5);
  std::cout << initGuess << "\n";

  auto container =
      random_vectors_from_guess<double, std::uniform_real_distribution>()(
          4, initGuess);
  std::cout << "Before sorting according to function value: \n";
  for (auto const &v : container) {
    std::cout << v << "|\n";
  }
  std::cout << "\n";
  std::cout << "After sorting according to function value: \n";

  auto f = [](vector_arg_t<double> const &args) {
    return (args(0) + args(1) + args(2));
  };

  std::sort(container.begin(), container.end(),
            [&](vector_arg_t<double> const &first,
                vector_arg_t<double> const &second) {
              return (f(first) < f(second));
            });

  for (auto const &v : container) {
    std::cout << v << "|\n";
  }
  std::cout << "\n";

  vector_arg_t<double> vat = vector_arg_t<double>::Zero(3);
  std::cout << vat << "\n";
}

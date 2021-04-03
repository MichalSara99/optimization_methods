#if !defined(OM_TEST_HELPERS)
#define OM_TEST_HELPERS

#include <iostream>

#include "utilities/om_types.hpp"

/**
 * @brief Contains test helpers
 *
 */

namespace om_test_helpers {

using om_types::f_vector_t;
using om_types::vector_arg_t;
using om_types::vector_t;
/**
 * @brief Helper for optimisation methods
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 */
template <typename fp_type = double> struct minimizer_helper {
  std::string name;
  f_vector_t<fp_type> objective;
  vector_arg_t<fp_type> guess;
  vector_arg_t<fp_type> minimizer;

  minimizer_helper(std::string const &name, f_vector_t<fp_type> objective,
                   vector_arg_t<fp_type> const &guess,
                   vector_arg_t<fp_type> const &minimizer)
      : objective(objective), guess(guess), name(name), minimizer(minimizer) {}
};
} // namespace om_test_helpers

#endif /// OM_TEST_HELPERS
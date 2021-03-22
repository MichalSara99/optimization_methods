#if !defined(OM_TEST_HELPERS)
#define OM_TEST_HELPERS

#include <iostream>

#include "utilities/om_types.hpp"

namespace om_test_helpers {

using om_types::f_vector_t;
using om_types::vector_arg_t;
using om_types::vector_t;

template <typename T = double> struct minimizer_helper {
  std::string name;
  f_vector_t<T> objective;
  vector_arg_t<T> guess;
  vector_arg_t<T> minimizer;

  minimizer_helper(std::string const &name, f_vector_t<T> objective,
                   vector_arg_t<T> const &guess,
                   vector_arg_t<T> const &minimizer)
      : objective(objective), guess(guess), name(name), minimizer(minimizer) {}
};
} // namespace om_test_helpers

#endif /// OM_TEST_HELPERS
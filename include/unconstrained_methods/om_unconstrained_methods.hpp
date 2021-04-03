#if !defined(OM_UNCONSTRAINED_METHODS)
#define OM_UNCONSTRAINED_METHODS

#include "multi_dim/conjugate_gradient/om_conjugate_gradient.hpp"
#include "multi_dim/quasi_newton/om_quasi_newton.hpp"
#include "multi_dim/steepest_descent/om_steepest_descent.hpp"
#include "multi_dim/zero_order/om_zero_order.hpp"
#include "one_dim/om_line_methods.hpp"
#include "utilities/om_types.hpp"
#include "utilities/om_utilities.hpp"

namespace om_unconstrained_methods {

using om_types::f_scalar_t;
using om_utilities::range;

/**
 * @brief Minimize one dimensional objective function
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 * @tparam line_search_method is any of the one_dim methods
 * @param objective objective function of f_scalar_t type
 * @param range range where to look for minimum
 * @param tolerance tolerance of minimiser
 * @param max_iters maximum number of iterations
 * @return std::tuple<fp_type, fp_type, std::size_t, std::size_t>
 */
template <typename fp_type = double,
          template <typename, typename>
          typename line_search_method = om_line_methods::brent_method>
std::tuple<fp_type, fp_type, std::size_t, std::size_t>
minimize(f_scalar_t<fp_type> &&objective, range<fp_type> const &range,
         fp_type tolerance, std::size_t const &max_iters) {
  line_search_method<fp_type, void> minimisor{range, tolerance, max_iters};
  return minimisor(std::forward<f_scalar_t<fp_type>>(objective));
}

} // namespace om_unconstrained_methods

#endif /// OM_UNCONSTRAINED_METHODS
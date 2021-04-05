#if !defined(OM_UNCONSTRAINED_METHODS)
#define OM_UNCONSTRAINED_METHODS

#include "multi_dim/conjugate_gradient/om_conjugate_gradient.hpp"
#include "multi_dim/quasi_newton/om_quasi_newton.hpp"
#include "multi_dim/steepest_descent/om_steepest_descent.hpp"
#include "multi_dim/zero_order/om_zero_order.hpp"
#include "om_unconstrained_methods_traits.hpp"
#include "one_dim/om_line_methods.hpp"
#include "utilities/om_types.hpp"
#include "utilities/om_utilities.hpp"

namespace om_unconstrained_methods {

using om_types::f_scalar_t;
using om_types::f_vector_t;
using om_types::vector_arg_t;
using om_types::vector_t;
using om_unconstrained_methods_traits::is_zero_order_method;
using om_utilities::range;

/**
 * @brief Minimize scalar objective function of one variable
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

/**
 * @brief Minimize scalar objective function of more then one variable (excluded
 * zero-order methods)
 * @details Zero-order methods (Nelder-Mead and Powell conjugate) are not
 * allowed here (this is taken care of via std::enable_if and
 * is_zero_order_method trait)
 *
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 * @tparam method optimisation method
 * @tparam line_search_method line search method
 * @tparam std::enable_if<
 * is_zero_order_method<method<fp_type>>::value>::type
 * @param objective objective function
 * @param init_guess initial guess
 * @param max_iters maximum number of iterations
 * @param arg_tol tolerance for a stopping criteria
 * @param grad_tol tolerance for gradient
 * @param fun_tol tolerance for a vlaue of function
 * @param line_search_range range for line search method
 * @return std::tuple<vector_t<fp_type>, fp_type, std::size_t>
 */
template <typename fp_type = double,
          template <typename> typename method =
              om_quasi_newton::broyden_fletcher_goldfarb_shanno_method,
          template <typename, typename>
          typename line_search_method = om_line_methods::golden_section_method,
          typename = typename std::enable_if<
              is_zero_order_method<method<fp_type>>::value>::type>
std::tuple<vector_t<fp_type>, fp_type, std::size_t>
minimize(f_vector_t<fp_type> &&objective,
         vector_arg_t<fp_type> const &init_guess, std::size_t const &max_iters,
         fp_type arg_tol = 1e-4, fp_type grad_tol = 1e-4,
         fp_type fun_tol = 1e-4,
         range<fp_type> const &line_search_range = range<fp_type>(-1.0, 1.0)) {

  line_search_method<fp_type, void> minimisor{line_search_range, arg_tol,
                                              max_iters};
  method<fp_type> m{minimisor, max_iters, arg_tol, grad_tol, fun_tol};
  return m.minimize(std::forward<f_vector_t<fp_type>>(objective), init_guess);
}

} // namespace om_unconstrained_methods

#endif /// OM_UNCONSTRAINED_METHODS
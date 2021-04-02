#if !defined(OM_FLETCHER_REEEVES)
#define OM_FLETCHER_REEEVES

#include "unconstrained_methods/multi_dim/conjugate_gradient/om_conjugate_gradient_base.hpp"
#include "unconstrained_methods/one_dim/om_line_methods.hpp"
#include "utilities/om_differentiation.hpp"
#include "utilities/om_types.hpp"

namespace om_unconstrained_methods {

namespace om_conjugate_gradient {

using om_differentiation::central_difference;
using om_types::f_vector_t;
using om_types::vector_arg_t;
using om_types::vector_t;
/**
 * @brief Fletcher-Reeves method object
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 */
template <typename fp_type = double>
class fletcher_reeves_method : public conjugate_gradient_base<fp_type> {
public:
  /**
   * @brief Construct a new fletcher reeves method object
   *
   * @param line_search_minimiser line method to be used in finding the
   * minimiser
   * @param max_iters maximum number of iterations
   * @param arg_tol tolerance for stopping criteria
   * @param grad_tol tolerance for gradient
   * @param fun_tol tolerance for a value of objective function
   */
  explicit fletcher_reeves_method(
      f_line_minimiser_t<fp_type> const &line_search_minimiser,
      std::size_t const &max_iters = 100, fp_type arg_tol = 1e-4,
      fp_type grad_tol = 1e-4, fp_type fun_tol = 1e-4)
      : conjugate_gradient_base<fp_type>{line_search_minimiser, max_iters,
                                         arg_tol, grad_tol, fun_tol} {}
  /**
   * @brief Function method that minimises the objective function
   *
   * @param objective objective function
   * @param init_guess initial guess
   * @return std::tuple<vector_t<fp_type>, fp_type, std::size_t>
   */
  std::tuple<vector_t<fp_type>, fp_type, std::size_t>
  minimize(f_vector_t<fp_type> objective,
           vector_arg_t<fp_type> const &init_guess) const;
};
} // namespace om_conjugate_gradient
} // namespace om_unconstrained_methods

template <typename fp_type>
std::tuple<om_unconstrained_methods::om_conjugate_gradient::vector_t<fp_type>,
           fp_type, std::size_t>
om_unconstrained_methods::om_conjugate_gradient::
    fletcher_reeves_method<fp_type>::minimize(
        om_unconstrained_methods::om_conjugate_gradient::f_vector_t<fp_type>
            objective,
        om_unconstrained_methods::om_conjugate_gradient::vector_arg_t<
            fp_type> const &init_guess) const {

  auto x_prev = init_guess;
  vector_t<fp_type> grad_prev;
  vector_t<fp_type> grad;
  vector_t<fp_type> u;
  vector_t<fp_type> x;
  fp_type beta{};
  std::size_t iters{0};

  // Start of Algorithm:
  grad_prev = central_difference<1, fp_type>()(objective, x_prev);
  u = static_cast<fp_type>(-1.0) * grad_prev;

  while (iters < this->max_iters_) {

    auto F_lambda = [&](fp_type const &lambda) {
      return objective(x_prev + lambda * u);
    };
    auto optim_lambda = this->lsm_(F_lambda);
    x = x_prev + (std::get<0>(optim_lambda)) * u;
    grad = central_difference<1, fp_type>()(objective, x);
    // Stopping criteria:
    if (((x - x_prev).norm() < this->arg_tol_) &&
        (grad.norm() < this->grad_tol_) &&
        (std::abs(objective(x) - objective(x_prev)) < this->fun_tol_)) {
      return std::make_tuple(x, objective(x), iters);
    } else {
      beta =
          ((grad.norm() / grad_prev.norm()) * (grad.norm() / grad_prev.norm()));
      u = (((static_cast<fp_type>(-1.0)) * grad) + (beta * u));
      x_prev = x;
      grad_prev = grad;
      iters++;
    }
  }
  return std::make_tuple(x, objective(x), iters);
}

#endif // OM_FLETCHER_REEEVES
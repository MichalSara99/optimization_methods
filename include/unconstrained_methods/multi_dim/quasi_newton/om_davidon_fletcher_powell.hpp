#if !defined(OM_DAVIDON_FLETCHER_POWELL)
#define OM_DAVIDON_FLETCHER_POWELL

#include "unconstrained_methods/multi_dim/quasi_newton/om_quasi_newton_base.hpp"
#include "unconstrained_methods/one_dim/om_line_methods.hpp"
#include "utilities/om_differentiation.hpp"
#include "utilities/om_types.hpp"

namespace om_unconstrained_methods {

namespace om_quasi_newton {

using om_differentiation::central_difference;
using om_types::f_vector_t;
using om_types::matrix_t;
using om_types::vector_arg_t;
using om_types::vector_t;

/**
 * @brief Davidon-Fletcher-Powell method object
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 */
template <typename fp_type = double>
class davidon_fletcher_powell_method : public quasi_newton_base<fp_type> {
public:
  /**
   * @brief Construct a new davidon fletcher powell method object
   *
   * @param line_search_minimiser line method to be used in finding the
   * minimiser
   * @param max_iters maximum number of iterations
   * @param arg_tol tolerance for stopping criteria
   * @param grad_tol tolerance for gradient
   * @param fun_tol tolerance for a value of objective function
   */
  davidon_fletcher_powell_method(
      f_line_minimiser_t<fp_type> const &line_search_minimiser,
      std::size_t const &max_iters = 100, fp_type arg_tol = 1e-4,
      fp_type grad_tol = 1e-4, fp_type fun_tol = 1e-4)
      : quasi_newton_base<fp_type>{line_search_minimiser, max_iters, arg_tol,
                                   grad_tol, fun_tol} {}
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
} // namespace om_quasi_newton
} // namespace om_unconstrained_methods

template <typename fp_type>
std::tuple<om_unconstrained_methods::om_quasi_newton::vector_t<fp_type>,
           fp_type, std::size_t>
om_unconstrained_methods::om_quasi_newton::
    davidon_fletcher_powell_method<fp_type>::minimize(
        om_unconstrained_methods::om_quasi_newton::f_vector_t<fp_type>
            objective,
        om_unconstrained_methods::om_quasi_newton::vector_arg_t<fp_type> const
            &init_guess) const {

  auto x_prev = init_guess;
  matrix_t<fp_type> A;
  matrix_t<fp_type> B;
  matrix_t<fp_type> G =
      Eigen::MatrixXd::Identity(init_guess.rows(), init_guess.rows());
  vector_t<fp_type> grad_prev;
  vector_t<fp_type> grad;
  vector_t<fp_type> u;
  vector_t<fp_type> x;
  vector_t<fp_type> v;
  vector_t<fp_type> y;
  std::size_t iters{0};

  // Start of Algorithm:
  grad_prev = central_difference<1, fp_type>()(objective, x_prev);
  u = static_cast<fp_type>(-1.0) * G * grad_prev;
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
      v = (std::get<0>(optim_lambda)) * u;
      y = grad - grad_prev;
      // rank 1-update:
      A = (v * v.transpose()) / (v.transpose() * y);
      B = static_cast<fp_type>(-1.0) * (G * y * (G * y).transpose()) /
          (y.transpose() * G * y);
      // update:
      G += (A + B);
      grad_prev = grad;
      x_prev = x;
      u = static_cast<fp_type>(-1.0) * G * grad_prev;
      iters++;
    }
  }
  return std::make_tuple(x, objective(x), iters);
}

#endif /// OM_DAVIDON_FLETCHER_POWELL
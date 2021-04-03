#if !defined(OM_STEEPEST_DESCENT)
#define OM_STEEPEST_DESCENT

#include "unconstrained_methods/one_dim/om_line_methods.hpp"
#include "utilities/om_differentiation.hpp"
#include "utilities/om_types.hpp"

namespace om_unconstrained_methods {
/**
 * @brief Contains steepest-descent method
 *
 */
namespace om_steepest_descent {

using om_differentiation::central_difference;
using om_types::f_line_minimiser_t;
using om_types::f_vector_t;
using om_types::vector_arg_t;
using om_types::vector_t;
/**
 * @brief Steepest descent method object
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 */
template <typename fp_type = double> class steepest_descent_method {
private:
  fp_type arg_tol_;
  fp_type grad_tol_;
  fp_type fun_tol_;
  std::size_t max_iters_;
  f_line_minimiser_t<fp_type> lsm_;

public:
  /**
   * @brief Construct a new steepest descent method object
   *
   * @param line_search_minimiser line method to be used in finding the
   * minimiser
   * @param max_iters maximum number of iterations
   * @param arg_tol tolerance for stopping criteria
   * @param grad_tol tolerance for gradient
   * @param fun_tol tolerance for a value of objective function
   */
  explicit steepest_descent_method(
      f_line_minimiser_t<fp_type> const &line_search_minimiser,
      std::size_t const &max_iters = 100, fp_type arg_tol = 1e-4,
      fp_type grad_tol = 1e-4, fp_type fun_tol = 1e-4)
      : arg_tol_{arg_tol}, fun_tol_{fun_tol}, grad_tol_{grad_tol},
        lsm_{line_search_minimiser}, max_iters_{max_iters} {}

  virtual ~steepest_descent_method() {}
  /**
   * @brief Copy constructor of a steepest descent method object
   *
   * @param copy copy is the object which we want to make a copy of
   */
  steepest_descent_method(steepest_descent_method const &copy)
      : arg_tol_{copy.arg_tol_}, fun_tol_{copy.fun_tol_},
        grad_tol_{copy.grad_tol_}, lsm_{copy.lsm_}, max_iters_{
                                                        copy.max_iters_} {}
  /**
   * @brief Assignment operator of a steepest descent method object
   *
   * @param copy
   * @return steepest_descent_method&
   */
  steepest_descent_method &operator=(steepest_descent_method const &copy) {
    if (&copy != this) {
      arg_tol_ = copy.arg_tol_;
      fun_tol_ = copy.fun_tol_;
      grad_tol_ = copy.grad_tol_;
      lsm_ = copy.lsm_;
      max_iters_ = copy.max_iters_;
    }
    return *this;
  }
  /**
   * @brief Set the stopping criteria tolerance object
   *
   * @param arg_tol tolerance for stopping criteria
   */
  inline void set_arg_tolerance(fp_type arg_tol) { arg_tol_ = arg_tol; }
  /**
   * @brief Set the fun tolerance object
   *
   * @param fun_tol tolerance for a value of function
   */
  inline void set_fun_tolerance(fp_type fun_tol) { fun_tol_ = fun_tol; }
  /**
   * @brief Set the grad tolerance object
   *
   * @param grad_tol tolerance for gradient
   */
  inline void set_grad_tolerance(fp_type grad_tol) { grad_tol_ = grad_tol; }
  /**
   * @brief Set the max iterations object
   *
   * @param iters maximum number of iterations
   */
  inline void set_max_iterations(std::size_t const &iters) {
    max_iters_ = iters;
  }
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

} // namespace om_steepest_descent
} // namespace om_unconstrained_methods

template <typename fp_type>
std::tuple<om_unconstrained_methods::om_steepest_descent::vector_t<fp_type>,
           fp_type, std::size_t>
om_unconstrained_methods::om_steepest_descent::steepest_descent_method<
    fp_type>::minimize(om_types::f_vector_t<fp_type> objective,
                       om_types::vector_arg_t<fp_type> const &init_guess)
    const {

  auto x_prev = init_guess;
  vector_t<fp_type> grad_prev;
  vector_t<fp_type> grad;
  vector_t<fp_type> u;
  vector_t<fp_type> x;
  std::size_t iters{0};

  grad_prev = central_difference<1, fp_type>()(objective, x_prev);
  u = grad_prev / (static_cast<fp_type>(-1.0) * grad_prev.norm());

  while (iters < max_iters_) {

    auto F_lambda = [&](fp_type const &lambda) {
      return objective(x_prev + lambda * u);
    };
    auto optim_lambda = lsm_(F_lambda);
    x = x_prev + (std::get<0>(optim_lambda)) * u;
    grad = central_difference<1, fp_type>()(objective, x);
    if (((x - x_prev).norm() < arg_tol_) || (grad.norm() < grad_tol_) ||
        (std::abs(objective(x) - objective(x_prev)) < fun_tol_)) {
      return std::make_tuple(x, objective(x), iters);
    } else {
      grad_prev = grad;
      x_prev = x;
      u = grad_prev / (static_cast<fp_type>(-1.0) * grad_prev.norm());
      iters++;
    }
  }
  return std::make_tuple(x, objective(x), iters);
}
#endif //_STEEPEST_DESCENT_METHODS
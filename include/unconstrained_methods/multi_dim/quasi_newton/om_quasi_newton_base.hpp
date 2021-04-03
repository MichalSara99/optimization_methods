#if !defined(OM_QUASI_NEWTON_BASE)
#define OM_QUASI_NEWTON_BASE

#include "utilities/om_types.hpp"

namespace om_unconstrained_methods {

namespace om_quasi_newton {

using om_types::f_line_minimiser_t;

/**
 * @brief Quasi-Newton base class
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 * @tparam std::enable_if<
 * std::is_floating_point<fp_type>::value>::type
 */
template <typename fp_type = double,
          typename = typename std::enable_if<
              std::is_floating_point<fp_type>::value>::type>
class quasi_newton_base {
protected:
  fp_type arg_tol_;
  fp_type grad_tol_;
  fp_type fun_tol_;
  std::size_t max_iters_;
  f_line_minimiser_t<fp_type> lsm_;

public:
  /**
   * @brief Construct a new quasi newton base object
   *
   * @param line_search_minimiser line method to be used in finding the
   * minimiser
   * @param max_iters maximum number of iterations
   * @param arg_tol tolerance for stopping criteria
   * @param grad_tol tolerance for gradient
   * @param fun_tol tolerance for a value of objective function
   */
  quasi_newton_base(f_line_minimiser_t<fp_type> const &line_search_minimiser,
                    std::size_t const &max_iters = 100, fp_type arg_tol = 1e-4,
                    fp_type grad_tol = 1e-4, fp_type fun_tol = 1e-4)
      : lsm_{line_search_minimiser}, max_iters_{max_iters}, arg_tol_{arg_tol},
        fun_tol_{fun_tol}, grad_tol_{grad_tol} {}

  virtual ~quasi_newton_base() {}
  /**
   * @brief Construct a new quasi newton base object
   *
   * @param copy copy is the object which we want to make a copy of
   */
  quasi_newton_base(quasi_newton_base const &copy)
      : lsm_{copy.lsm_}, max_iters_{copy.max_iters_}, arg_tol_{copy.arg_tol_},
        fun_tol_{copy.fun_tol_}, grad_tol_{copy.grad_tol_} {}
  /**
   * @brief Assignment operator of a quasi newton base object
   *
   * @param copy
   * @return quasi_newton_base&
   */
  quasi_newton_base &operator=(quasi_newton_base const &copy) {
    if (this != &copy) {
      lsm_ = copy.lsm_;
      arg_tol_ = copy.arg_tol_;
      fun_tol_ = copy.fun_tol_;
      grad_tol_ = copy.grad_tol_;
      max_iters_ = copy.max_iters_;
    }
    return *this;
  }
  /**
   * @brief Set the max iterations object
   *
   * @param iters maximum number of iterations
   */
  inline void set_max_iterations(std::size_t const &iters) {
    max_iters_ = iters;
  }
  /**
   * @brief  Set the stopping criteria tolerance object
   *
   * @param arg_tol tolerance for stopping criteria
   */
  inline void set_arg_tolerance(fp_type arg_tol) { arg_tol_ = arg_tol; }
  /**
   * @brief Set the fun tolerance object
   *
   * @param fun_tol tolerance for a value of objective function
   */
  inline void set_fun_tolerance(fp_type fun_tol) { fun_tol_ = fun_tol; }
  /**
   * @brief Set the grad tolerance object
   *
   * @param grad_tol tolerance for gardient
   */
  inline void set_grad_tolerance(fp_type grad_tol) { grad_tol_ = grad_tol; }
};

} // namespace om_quasi_newton

} // namespace om_unconstrained_methods

#endif /// OM_QUASI_NEWTON_BASE
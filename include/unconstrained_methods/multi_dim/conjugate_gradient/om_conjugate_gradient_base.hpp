#if !defined(OM_CONJUGATE_GRADIENT_BASE)
#define OM_CONJUGATE_GRADIENT_BASE

#include "utilities/om_types.hpp"

namespace om_unconstrained_methods {

namespace om_conjugate_gradient {

using om_types::f_line_minimiser_t;

template <typename fp_type = double,
          typename = typename std::enable_if<
              std::is_floating_point<fp_type>::value>::type>
class conjugate_gradient_base {
protected:
  fp_type arg_tol_;
  fp_type grad_tol_;
  fp_type fun_tol_;
  std::size_t max_iters_;
  f_line_minimiser_t<fp_type> lsm_;

public:
  explicit conjugate_gradient_base(
      f_line_minimiser_t<fp_type> const &line_search_minimiser,
      std::size_t const &max_iters = 100, fp_type arg_tol = 1e-4,
      fp_type grad_tol = 1e-4, fp_type fun_tol = 1e-4)
      : arg_tol_{arg_tol}, lsm_{line_search_minimiser},
        max_iters_{max_iters}, grad_tol_{grad_tol}, fun_tol_{fun_tol} {}

  virtual ~conjugate_gradient_base() {}

  conjugate_gradient_base(conjugate_gradient_base const &copy)
      : arg_tol_{copy.arg_tol_}, lsm_{copy.lsm_}, max_iters_{copy.max_iters_},
        grad_tol_{copy.grad_tol_}, fun_tol_{copy.fun_tol_} {}

  conjugate_gradient_base &operator=(conjugate_gradient_base const &copy) {
    if (this != &copy) {
      arg_tol_ = copy.arg_tol_;
      fun_tol_ = copy.fun_tol_;
      grad_tol_ = copy.grad_tol_;
      lsm_ = copy.lsm_;
      max_iters_ = copy.max_iters_;
    }
    return *this;
  }

  inline void set_arg_tolerance(fp_type arg_tol) { arg_tol_ = arg_tol; }
  inline void set_fun_tolerance(fp_type fun_tol) { fun_tol_ = fun_tol; }
  inline void set_grad_tolerance(fp_type grad_tol) { grad_tol_ = grad_tol; }
  inline void set_max_iterations(std::size_t const &iters) {
    max_iters_ = iters;
  }
};

} // namespace om_conjugate_gradient
} // namespace om_unconstrained_methods
#endif // OM_CONJUGATE_GRADIENT_BASE
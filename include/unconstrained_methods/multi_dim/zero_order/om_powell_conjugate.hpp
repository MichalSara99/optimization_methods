#if !defined(OM_POWELL_CONJUGATE)
#define OM_POWELL_CONJUGATE

#include "unconstrained_methods/one_dim/om_line_methods.hpp"
#include "utilities/om_differentiation.hpp"
#include "utilities/om_types.hpp"

namespace om_unconstrained_methods {

namespace om_zero_order {

using om_types::f_line_minimiser_t;
using om_types::f_vector_t;
using om_types::vector_arg_t;
using om_types::vector_t;
using om_utilities::cartesian_basis_vectors;

/**
 * @brief Powell conjugate method object
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 */
template <typename fp_type = double> class powell_conjugate_method {
private:
  fp_type converge_tol_;
  std::size_t max_iters_;
  f_line_minimiser_t<fp_type> lsm_;

public:
  /**
   * @brief Construct a new powell conjugate method object
   *
   * @param line_search_minimiser line method to be used in finding the
   * minimiser
   * @param max_iters maximum number of iterations
   * @param convergence_tol tolerance for convergance
   */
  powell_conjugate_method(
      f_line_minimiser_t<fp_type> const &line_search_minimiser,
      std::size_t const &max_iters = 50, fp_type convergence_tol = 10e-4)
      : lsm_{line_search_minimiser}, max_iters_{max_iters},
        converge_tol_{convergence_tol} {}

  virtual ~powell_conjugate_method() {}
  /**
   * @brief Copy constructor a new powell conjugate method object
   *
   * @param copy copy is the object which we want to make a copy of
   */
  powell_conjugate_method(powell_conjugate_method const &copy)
      : max_iters_{copy.max_iters_}, lsm_{copy.lsm_}, converge_tol_{
                                                          copy.converge_tol_} {}

  /**
   * @brief Assignment operator of a powell conjugate method object
   *
   * @param copy
   * @return powell_conjugate_method&
   */
  powell_conjugate_method &operator=(powell_conjugate_method const &copy) {
    if (this != &copy) {
      lsm_ = copy.lsm_;
      max_iters_ = copy.max_iters_;
      converge_tol_ = copy.converge_tol_;
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
   * @brief Set the converge tolerance object
   *
   * @param converge_tol tolerance for convergance
   */
  inline void set_converge_tolerance(double converge_tol) {
    converge_tol_ = converge_tol;
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
} // namespace om_zero_order

template <typename fp_type>
std::tuple<om_zero_order::vector_t<fp_type>, fp_type, std::size_t>
om_zero_order::powell_conjugate_method<fp_type>::minimize(
    om_zero_order::f_vector_t<fp_type> objective,
    om_zero_order::vector_arg_t<fp_type> const &init_guess) const {

  auto x_prev = init_guess;
  std::size_t dim = init_guess.rows();
  std::size_t iters{0};

  auto u = cartesian_basis_vectors<fp_type>()(dim);
  auto x = std::vector<vector_t<fp_type>>(dim + 2);
  x[0] = x_prev;
  x[1] = x_prev;

  while (iters < this->max_iters_) {
    for (auto t = 0; t < u.size(); ++t) {

      auto F_lambda = [&](fp_type const &lambda) {
        return objective(x[t + 1] + lambda * u[t]);
      };
      auto optim_lambda = this->lsm_(F_lambda);
      x[t + 2] = x[t + 1] + std::get<0>(optim_lambda) * u[t];
    }

    for (auto t = 0; t < u.size() - 1; ++t) {
      u[t] = u[t + 1];
    }
    u[u.size() - 1] = x.back() - x[1];

    if ((x.back() - x[1]).norm() <= this->converge_tol_) {
      return std::make_tuple(x.back(), objective(x.back()), iters);
    } else {
      iters++;
      x[0] = x.back();
    }

    if ((iters % (dim + 1)) == 0) {
      u = cartesian_basis_vectors<fp_type>()(dim);
      x[1] = x.front();
    }
  }
  return std::make_tuple(x.back(), objective(x.back()), iters);
}

} // namespace om_unconstrained_methods
#endif /// OM_POWELL_CONJUGATE
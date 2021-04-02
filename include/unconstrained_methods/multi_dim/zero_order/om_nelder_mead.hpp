#if !defined(OM_NELDER_MEAD)
#define OM_NELDER_MEAD

#include "unconstrained_methods/one_dim/om_line_methods.hpp"
#include "utilities/om_differentiation.hpp"
#include "utilities/om_types.hpp"

namespace om_unconstrained_methods {

namespace om_nelder_mead {

using om_types::f_vector_t;
using om_types::vector_arg_t;
using om_types::vector_t;
using om_utilities::random_vectors_from_guess;

template <typename fp_type = double> class nelder_mead_method {
private:
  fp_type reflection_rho_;
  fp_type expansion_rho_;
  fp_type contraction_rho_;
  fp_type shrinkage_rho_;
  std::size_t max_iters_;
  fp_type converge_tol_;

  void check_rhos();

public:
  nelder_mead_method(std::size_t const &max_iters = 80,
                     fp_type convergence_tol = 10e-4,
                     fp_type reflection_rho = 0.5, fp_type expansion_rho = 1.5,
                     fp_type contraction_rho = 0.25,
                     fp_type shrinkage_rho = 0.5)
      : max_iters_{max_iters}, converge_tol_{convergence_tol},
        reflection_rho_{reflection_rho}, expansion_rho_{expansion_rho},
        contraction_rho_{contraction_rho}, shrinkage_rho_{shrinkage_rho} {
    check_rhos();
  }

  virtual ~nelder_mead_method() {}

  nelder_mead_method(nelder_mead_method const &copy)
      : max_iters_{copy.max_iters_}, converge_tol_{copy.converge_tol_},
        reflection_rho_{copy.reflection_rho_},
        expansion_rho_{copy.expansion_rho_},
        contraction_rho_{copy.contraction_rho_}, shrinkage_rho_{
                                                     copy.shrinkage_rho_} {}

  nelder_mead_method &operator=(nelder_mead_method const &copy) {
    if (this != &copy) {
      max_iters_ = copy.max_iters_;
      converge_tol_ = copy.converge_tol_;
      reflection_rho_ = copy.reflection_rho_;
      expansion_rho_ = copy.expansion_rho_;
      contraction_rho_ = copy.contraction_rho_;
      shrinkage_rho_ = copy.shrinkage_rho_;
    }
    return *this;
  }

  inline void set_max_iterations(std::size_t const &iters) {
    max_iters_ = iters;
  }
  inline void set_converge_tolerance(fp_type converge_tol) {
    converge_tol_ = converge_tol;
  }
  inline void set_reflection_rho(fp_type value) {
    assert((value > 0.0));
    reflection_rho_ = value;
  }
  inline void set_expansion_rho(fp_type value) {
    assert((value > 1.0));
    expansion_rho_ = value;
  }
  inline void set_contraction_rho(fp_type value) {
    assert(((value >= 0.0) && (value <= 0.5)));
    contraction_rho_ = value;
  }
  inline void set_shrinkage_rho(fp_type value) {
    assert(((value >= 0.0) && (value <= 1.0)));
    shrinkage_rho_ = value;
  }

  std::tuple<vector_t<fp_type>, fp_type, std::size_t>
  minimize(f_vector_t<fp_type> objective,
           vector_arg_t<fp_type> const &init_guess) const;
};

} // namespace om_nelder_mead

template <typename fp_type>
void om_nelder_mead::nelder_mead_method<fp_type>::check_rhos() {
  assert((this->reflection_rho_ > 0.0));
  assert((this->expansion_rho_ > 1.0));
  assert(((this->contraction_rho_ >= 0.0) && (this->contraction_rho_ <= 0.5)));
  assert(((this->shrinkage_rho_ >= 0.0) && (this->shrinkage_rho_ <= 1.0)));
}

template <typename fp_type>
std::tuple<vector_t<fp_type>, fp_type, std::size_t>
om_nelder_mead::nelder_mead_method<fp_type>::minimize(
    f_vector_t<fp_type> objective,
    vector_arg_t<fp_type> const &init_guess) const {

  std::size_t N = init_guess.rows();
  std::size_t iters{0};
  auto init_simplex =
      random_vectors_from_guess<fp_type, std::uniform_real_distribution>()(
          N + 1, init_guess);

  while (iters < this->max_iters_) {

    // 1.STEP: Sort:
    std::sort(init_simplex.begin(), init_simplex.end(),
              [&](vector_arg_t<fp_type> const &first,
                  vector_arg_t<fp_type> const &second) {
                return (objective(first) < objective(second));
              });
    // Convergence check:
    if ((init_simplex.front() - init_simplex.back()).norm() <=
        this->converge_tol_) {
      return std::make_tuple(init_simplex.front(),
                             objective(init_simplex.front()), iters);
    }
    // 2.STEP: Centroid:
    vector_arg_t<fp_type> centroid = vector_arg_t<fp_type>::Zero(N);
    for (auto t = 0; t < init_simplex.size() - 1; ++t) {
      centroid += ((1.0 / N) * init_simplex[t]);
    }
    // 3.STEP: Reflection:
    vector_arg_t<fp_type> reflection = vector_arg_t<fp_type>::Zero(N);
    reflection =
        centroid + this->reflection_rho_ * (centroid - init_simplex.back());

    if ((objective(centroid) <= objective(reflection)) &&
        (objective(reflection) <
         objective(init_simplex[init_simplex.size() - 2]))) {
      init_simplex[init_simplex.size() - 1] = reflection;
      iters++;
      // go to 1.STEP
      continue;
    } else if ((objective(reflection) < objective(centroid))) {
      // 4.STEP: Expansion:
      vector_arg_t<fp_type> expansion = vector_arg_t<fp_type>::Zero(N);
      expansion =
          centroid + this->expansion_rho_ * (centroid - init_simplex.back());

      if (objective(expansion) < objective(reflection)) {
        init_simplex[init_simplex.size() - 1] = expansion;
      } else {
        init_simplex[init_simplex.size() - 1] = reflection;
      }
      iters++;
      continue;
    } else {
      // 5.STEP: Contraction:
      vector_arg_t<fp_type> contraction = vector_arg_t<fp_type>::Zero(N);
      contraction =
          centroid + this->contraction_rho_ * (centroid - init_simplex.back());

      if (objective(contraction) < objective(init_simplex.back())) {
        init_simplex[init_simplex.size() - 1] = contraction;
        iters++;
        continue;
      } else {
        // 6.STEP: Shrinking:
        for (auto t = 1; t < init_simplex.size(); ++t) {
          init_simplex[t] =
              (init_simplex.front() +
               this->shrinkage_rho_ * (init_simplex[t] - init_simplex.front()));
        }
        iters++;
      }
    }
  }
  return std::make_tuple(init_simplex.front(), objective(init_simplex.front()),
                         iters);
}

} // namespace om_unconstrained_methods

#endif /// OM_NELDER_MEAD
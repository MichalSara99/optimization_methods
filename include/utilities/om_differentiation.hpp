#if !defined(OM_DIFFERENTIATION)
#define OM_DIFFERENTIATION

#include "om_differentiation_traits.hpp"
#include "om_types.hpp"

/**
 * @brief Contains some numerical differentiation functors
 *
 */
namespace om_differentiation {

// loding types:
using om_differentiation_traits::central_difference_trait;
using om_differentiation_traits::forward_difference_trait;
using om_types::f_scalar_t;
using om_types::f_vector_t;
using om_types::vector_arg_t;
using om_types::vector_t;

/**
 * @brief Divided difference functor
 *
 * @tparam order order of difference
 * @tparam fp_type fp_type is a floating-point template parameter
 * @details order = 0, order = 1, order = 2, order = 3 currently supported
 */
template <std::size_t order, typename fp_type = double,
          typename = typename std::enable_if<order >= 0 && order <= 3>::type>
struct divided_difference {};

template <typename fp_type> struct divided_difference<0, fp_type> {
  fp_type operator()(f_scalar_t<fp_type> fun, fp_type const &arg) const {
    return fun(arg);
  }
};

template <typename fp_type> struct divided_difference<1, fp_type> {
  fp_type operator()(f_scalar_t<fp_type> fun,
                     std::tuple<fp_type, fp_type> const &arg) const {
    return ((divided_difference<0, fp_type>()(fun, std::get<1>(arg)) -
             divided_difference<0, fp_type>()(fun, std::get<0>(arg))) /
            (std::get<1>(arg) - std::get<0>(arg)));
  }
};

template <typename fp_type> struct divided_difference<2, fp_type> {
  fp_type operator()(f_scalar_t<fp_type> fun,
                     std::tuple<fp_type, fp_type, fp_type> const &arg) const {
    return ((divided_difference<1, fp_type>()(
                 fun, std::make_tuple(std::get<1>(arg), std::get<2>(arg))) -
             divided_difference<1, fp_type>()(
                 fun, std::make_tuple(std::get<0>(arg), std::get<1>(arg)))) /
            (std::get<2>(arg) - std::get<0>(arg)));
  }
};

/**
 * @brief forward difference functor
 *
 * @tparam order order of difference
 * @tparam fp_type
 * @tparam std::enable_if<
 * std::is_floating_point<fp_type>::value>::type
 * @details order = 0,order = 1 currently supported
 */
template <std::size_t order, typename fp_type,
          typename = typename std::enable_if<
              std::is_floating_point<fp_type>::value>::type>
struct forward_difference;

template <typename fp_type> struct forward_difference<0, fp_type> {
public:
  vector_t<fp_type> operator()(f_vector_t<fp_type> fun,
                               vector_arg_t<fp_type> const &args) const {
    auto result = args;
    for (auto t = 0; t < args.rows(); ++t) {
      result(t) = fun(args);
    }
    return result;
  }
};

template <typename fp_type> struct forward_difference<1, fp_type> {

  vector_t<fp_type> operator()(f_vector_t<fp_type> fun,
                               vector_arg_t<fp_type> const &args) const {

    const fp_type step = forward_difference_trait<fp_type>::step_size;
    fp_type base = fun(args);
    auto result = args;
    for (auto t = 0; t < args.rows(); ++t) {
      auto a = args;
      a(t) += step;
      result(t) = (fun(a) - base) / step;
    }
    return result;
  }
};

/**
 * @brief central difference functor
 *
 * @tparam order order of difference
 * @tparam fp_type fp_type is a floating-point template parameter
 * @tparam std::enable_if<
 * std::is_floating_point<fp_type>::value>::type
 * @details order = 0,order = 1 currently supported
 */
template <std::size_t order, typename fp_type,
          typename = typename std::enable_if<
              std::is_floating_point<fp_type>::value>::type>
struct central_difference {};

template <typename fp_type> struct central_difference<0, fp_type> {
public:
  vector_t<fp_type> operator()(f_vector_t<fp_type> fun,
                               vector_arg_t<fp_type> const &args) const {
    auto result = args;
    for (auto t = 0; t < args.rows(); ++t) {
      result(t) = fun(args);
    }
    return result;
  }
};

template <typename fp_type> struct central_difference<1, fp_type> {
  vector_t<fp_type> operator()(f_vector_t<fp_type> fun,
                               vector_arg_t<fp_type> const &args) const {
    const fp_type step = central_difference_trait<fp_type>::step_size;
    auto result = args;
    for (auto t = 0; t < args.rows(); ++t) {
      auto a_plus = args;
      auto a_minus = args;
      a_minus(t) -= step;
      a_plus(t) += step;
      result(t) =
          (fun(a_plus) - fun(a_minus)) / (static_cast<fp_type>(2.0) * step);
    }
    return result;
  }
};

} // namespace om_differentiation

#endif // OM_DIFFERENTIATION
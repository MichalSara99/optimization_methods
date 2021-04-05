#if !defined(OM_UNCONSTRAINED_METHODS_TRAITS)
#define OM_UNCONSTRAINED_METHODS_TRAITS

#include "multi_dim/zero_order/om_zero_order.hpp"

/**
 * @brief Contains some traits for minimize() global function
 *
 */
namespace om_unconstrained_methods_traits {

template <typename method> struct is_zero_order_method {
  static const bool value = true;
};

template <>
struct is_zero_order_method<
    om_unconstrained_methods::om_zero_order::nelder_mead_method<>> {
  static const bool value = false;
};

template <>
struct is_zero_order_method<
    om_unconstrained_methods::om_zero_order::powell_conjugate_method<>> {
  static const bool value = false;
};

} // namespace om_unconstrained_methods_traits

#endif /// OM_UNCONSTRAINED_METHODS_TRAITS
#if !defined(OM_DIFFERENTIATION_TRAITS)
#define OM_DIFFERENTIATION_TRAITS

/**
 * @brief Contains traits tested for numerical differentiation
 *
 */
namespace om_differentiation_traits {

/**
 * @brief forward difference trait
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 */
template <typename fp_type> struct forward_difference_trait {
  static constexpr fp_type step_size = 10e-6;
};
/**
 * @brief central difference trait
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 */
template <typename fp_type> struct central_difference_trait {
  static constexpr fp_type step_size = 10e-7;
};

} // namespace om_differentiation_traits

#endif // OM_DIFFERENTIATION_TRAITS
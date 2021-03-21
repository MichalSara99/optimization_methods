#if !defined(OM_DIFFERENTIATION_TRAITS)
#define OM_DIFFERENTIATION_TRAITS

namespace om_differentiation_traits {

template <typename T> struct forward_difference_trait {
  static constexpr T step_size = 10e-6;
};

template <typename T> struct central_difference_trait {
  static constexpr T step_size = 10e-7;
};

} // namespace om_differentiation_traits

#endif // OM_DIFFERENTIATION_TRAITS
#if !defined(OM_COMMON)
#define OM_COMMON

#include "om_types.hpp"
#include <type_traits>
/**
 * @brief Contains some commonly used measures
 *
 */
namespace om_common {

using om_types::f_scalar_t;

/**
 * @brief furthest_from functor
 *
 * @tparam count number of points to measure the distance from
 * @tparam fp_type fp_type is a floating-point template parameter
 * @details currently count = 2, count = 3 is supported
 */
template <std::size_t count, typename fp_type = double,
          typename = typename std::enable_if<count >= 2 && count <= 3>::type>
struct furthest_from {};

template <typename fp_type> struct furthest_from<2, fp_type> {
public:
  fp_type operator()(fp_type const &target, fp_type const &x1,
                     fp_type const &x2) const {
    if (std::abs(target - x1) > std::abs(target - x2)) {
      return x1;
    } else {
      return x2;
    }
  }
};

template <typename fp_type> struct furthest_from<3, fp_type> {
public:
  fp_type operator()(fp_type const &target, fp_type const &x1,
                     fp_type const &x2, fp_type const &x3) const {
    auto first = furthest_from<2, fp_type>()(target, x1, x2);
    auto second = furthest_from<2, fp_type>()(target, x2, x3);
    if (std::abs(first - target) > std::abs(second - target)) {
      return first;
    } else {
      return second;
    }
  }
};

/**
 * @brief closest_to functor
 *
 * @tparam count number of points
 * @tparam fp_type fp_type is a floating-point template parameter
 * @details count = 2,count = 3 is currently supported
 */

template <std::size_t count, typename fp_type = double,
          typename = typename std::enable_if<count >= 2 && count <= 3>::type>
struct closest_to {};

template <typename fp_type> struct closest_to<2, fp_type> {
public:
  fp_type operator()(fp_type const &target, fp_type const &x1,
                     fp_type const &x2) const {
    if (std::abs(target - x1) < std::abs(target - x2)) {
      return x1;
    } else {
      return x2;
    }
  }
};

template <typename fp_type> struct closest_to<3, fp_type> {
public:
  fp_type operator()(fp_type const &target, fp_type const &x1,
                     fp_type const &x2, fp_type const &x3) const {
    auto first = closest_to<2, fp_type>()(target, x1, x2);
    auto second = closest_to<2, fp_type>()(target, x2, x3);
    if (std::abs(first - target) < std::abs(second - target)) {
      return first;
    } else {
      return second;
    }
  }
};

/**
 * @brief min_arg functor retuns argument as which a function takes minimum
 * value
 *
 * @tparam count number of arguments
 * @tparam fp_type fp_type is a floating-point template parameter
 * @details count = 2,count = 3 currently supported
 */
template <std::size_t count, typename fp_type = double,
          typename = typename std::enable_if<count >= 2 && count <= 3>::type>
struct min_arg {};

template <typename fp_type> struct min_arg<2, fp_type> {
public:
  std::pair<fp_type, fp_type> operator()(f_scalar_t<fp_type> fun,
                                         fp_type const &first,
                                         fp_type const &second) const {
    auto f_first = fun(first);
    auto f_second = fun(second);
    return (f_first < f_second ? std::make_pair(first, f_first)
                               : std::make_pair(second, f_second));
  }
};

template <typename fp_type> struct min_arg<3, fp_type> {
public:
  std::pair<fp_type, fp_type> operator()(f_scalar_t<fp_type> fun,
                                         fp_type const &first,
                                         fp_type const &second,
                                         fp_type const &third) const {
    auto tmp1 = min_arg<2, fp_type>()(fun, first, second);
    auto tmp2 = min_arg<2, fp_type>()(fun, second, third);
    return min_arg<2, fp_type>()(fun, tmp1.first, tmp2.first);
  }
};

/**
 * @brief max_arg functor returns argument at which a function takes maximum
 * value
 *
 * @tparam count number of arguments
 * @tparam fp_type fp_type is a floating-point template argument
 * @details count = 2,count = 3 currently supported
 */
template <std::size_t count, typename fp_type = double,
          typename = typename std::enable_if<count >= 2 && count <= 3>::type>
struct max_arg {};

template <typename fp_type> struct max_arg<2, fp_type> {
public:
  std::pair<fp_type, fp_type> operator()(f_scalar_t<fp_type> fun,
                                         fp_type const &first,
                                         fp_type const &second) const {
    fp_type f_first = fun(first);
    fp_type f_second = fun(second);
    return (f_first > f_second ? std::make_pair(first, f_first)
                               : std::make_pair(second, f_second));
  }
};

template <typename fp_type> struct max_arg<3, fp_type> {
public:
  std::pair<fp_type, fp_type> operator()(f_scalar_t<fp_type> fun,
                                         fp_type const &first,
                                         fp_type const &second,
                                         fp_type const &third) const {
    std::pair<fp_type, fp_type> tmp1 =
        max_arg<2, fp_type>()(fun, first, second);
    std::pair<fp_type, fp_type> tmp2 =
        max_arg<2, fp_type>()(fun, second, third);
    return max_arg<2, fp_type>()(fun, tmp1.first, tmp2.first);
  }
};

} // namespace om_common

#endif /// OM_COMMON
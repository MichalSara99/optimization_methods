#if !defined(OM_TEST_FUNCTIONS)
#define OM_TEST_FUNCTIONS

#include "om_test_helpers.hpp"
#include "utilities/om_macros.hpp"
#include "utilities/om_types.hpp"
#include <cmath>
#include <memory>

using namespace om_types;
using namespace om_test_helpers;
using namespace std::string_literals;

/**
 * @brief Some classical test functions (designed by Rao)
 *
 */
namespace om_test_functions {

/**
 * @brief Pi definition used in the Rao test functons
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 */
template <typename fp_type> constexpr fp_type pi{3.14159265359};

/**
 * @brief Rosenbrock's parabolic valley test function
 * @details initial guess = (-1.2,1.0),
 * minimiser = (1.0,1.0)
 *
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 * @param args arguments of the function
 * @return fp_type
 */
template <typename fp_type>
fp_type rosenbrock_parabolic_valley(vector_arg_t<fp_type> const &args) {
  return (static_cast<fp_type>(100.0) * (args(1) - args(0) * args(0)) *
              (args(1) - args(0) * args(0)) +
          (static_cast<fp_type>(1.0) - args(0)) *
              (static_cast<fp_type>(1.0) - args(0)));
}

/**
 * @brief Quadratic test function
 * @details initial guess = (0.0,0.0), minimiser = (1.0,3.0)
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 * @param args function arguments
 * @return fp_type
 */
template <typename fp_type>
fp_type quadratic_function(vector_arg_t<fp_type> const &args) {
  return ((args(0) + static_cast<fp_type>(2.0) * args(1) -
           static_cast<fp_type>(7.0)) *
              (args(0) + static_cast<fp_type>(2.0) * args(1) -
               static_cast<fp_type>(7.0)) +
          (static_cast<fp_type>(2.0) * args(0) + args(1) -
           static_cast<fp_type>(5.0)) *
              (static_cast<fp_type>(2.0) * args(0) + args(1) -
               static_cast<fp_type>(5.0)));
}

/**
 * @brief Powell's quadratic test function
 * @details initial guess = (3.0,-1.0,0.0,1.0),
 * minimiser = (0.0,0.0,0.0,0.0)
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 * @param args function arguments
 * @return fp_type
 */
template <typename fp_type>
fp_type powell_function(vector_arg_t<fp_type> const &args) {
  return ((args(0) + static_cast<fp_type>(10.0) * args(1)) *
              (args(0) + static_cast<fp_type>(10.0) * args(1)) +
          static_cast<fp_type>(5.0) * (args(2) - args(3)) *
              (args(2) - args(3)) +
          std::pow<fp_type>((args(1) - static_cast<fp_type>(2.0) * args(2)),
                            static_cast<fp_type>(4.0)) +
          std::pow<fp_type>(static_cast<fp_type>(10.0) * (args(0) - args(3)),
                            static_cast<fp_type>(4.0)));
}

/**
 * @brief Fletcher and Powell's helical valley test function
 * @details initial guess = (-1.0,0.0,0.0),
 * minimiser = (1.0,0.0,0.0)
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 * @param args function arguments
 * @return fp_type
 */
template <typename fp_type>
fp_type fletcher_powell_helical_valley(vector_arg_t<fp_type> const &args) {
  fp_type theta{0.0};
  if (args(0) > static_cast<fp_type>(0.0)) {
    theta =
        (static_cast<fp_type>(1.0) / (static_cast<fp_type>(2.0) * pi<fp_type>)*(
                                         std::atan(args(1) / args(0))));
  } else {
    theta = (static_cast<fp_type>(1.0) /
             (static_cast<fp_type>(2.0) *
              pi<fp_type>)*(pi<fp_type> + (std::atan(args(1) / args(0)))));
  }
  return (static_cast<fp_type>(100.0) *
              (((args(2) - static_cast<fp_type>(10.0) * theta) *
                (args(2) - static_cast<fp_type>(10.0) * theta)) +
               ((std::sqrt(args(0) * args(0) + args(1) * args(1)) -
                 static_cast<fp_type>(1.0)) *
                (std::sqrt(args(0) * args(0) + args(1) * args(1)) -
                 static_cast<fp_type>(1.0)))) +
          args(2) * args(2));
}

/**
 * @brief Non-linear test function of 3 variables
 * @details initial guess = (0.0,1.0,2.0),
 * minimiser = (1.0,1.0,1.0)
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 * @param args function arguments
 * @return fp_type
 */
template <typename fp_type>
fp_type non_linear_function(vector_arg_t<fp_type> const &args) {
  auto term_1 =
      static_cast<fp_type>(1.0) /
      (static_cast<fp_type>(1.0) + ((args(0) - args(1)) * (args(0) - args(1))));
  auto term_2 =
      std::sin(static_cast<fp_type>(0.5) * pi<fp_type> * args(1) * args(2));
  auto term_3 =
      std::exp(static_cast<fp_type>(-1.0) *
               (((args(0) + args(2)) / args(1)) - static_cast<fp_type>(2.0)) *
               (((args(0) + args(2)) / args(1)) - static_cast<fp_type>(2.0)));
  return (static_cast<fp_type>(-1.0) * term_1 - term_2 - term_3);
}

/**
 * @brief Freudenstein and Roth test function
 * @details initial guess = (0.5,-2.0),
 * minimiser = (5.0,4.0),
 * local_minimiser = (11.41..., -0.8968)
 *
 * @todo Check if the minimiser and local_minimiser are correct!!
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 * @param args function arguments
 * @return fp_type
 */
template <typename fp_type>
fp_type freudenstein_roth_function(vector_arg_t<fp_type> const &args) {
  auto term_1 = std::pow((static_cast<fp_type>(-13.0) + args(0) +
                          ((static_cast<fp_type>(5.0) - args(1)) * args(1) -
                           static_cast<fp_type>(2.0)) *
                              args(1)),
                         static_cast<fp_type>(2.0));
  auto term_2 = std::pow((static_cast<fp_type>(-29.0) + args(0) +
                          ((args(1) + static_cast<fp_type>(1.0)) * args(1) -
                           static_cast<fp_type>(14.0)) *
                              args(1)),
                         static_cast<fp_type>(2.0));
  return (term_1 + term_2);
}

/**
 * @brief Powell's badly scaled test function
 * @details initial guess = (0.0,1.0),
 * minimiser = (1.098...*10^-5,9.106...)
 * @todo Check the validity of minimiser!!!
 *
 * @tparam fp_type
 * @param args
 * @return fp_type
 */
template <typename fp_type>
fp_type powell_badly_scaled_function(vector_arg_t<fp_type> const &args) {
  auto term_1 = (static_cast<fp_type>(10'000) * args(0) * args(1) -
                 static_cast<fp_type>(1.0)) *
                (static_cast<fp_type>(10'000) * args(0) * args(1) -
                 static_cast<fp_type>(1.0));
  auto term_2 = (std::exp(static_cast<fp_type>(-1.0) * args(0)) +
                 std::exp(static_cast<fp_type>(-1.0) * args(1)) -
                 static_cast<fp_type>(1.0001)) *
                (std::exp(static_cast<fp_type>(-1.0) * args(0)) +
                 std::exp(static_cast<fp_type>(-1.0) * args(1)) -
                 static_cast<fp_type>(1.0001));
  return (term_1 + term_2);
}

/**
 * @brief Beale's test function
 * @details initial guess = (1.0,1.0),
 * minimiser = (3.0,0.5)
 *
 * @tparam fp_type
 * @param args
 * @return fp_type
 */
template <typename fp_type>
fp_type beale_function(vector_arg_t<fp_type> const &args) {
  auto term_1 = std::pow((static_cast<fp_type>(1.5) -
                          args(0) * (static_cast<fp_type>(1.0) - args(1))),
                         static_cast<fp_type>(2.0));
  auto term_2 =
      std::pow((static_cast<fp_type>(2.25) -
                args(0) * (static_cast<fp_type>(1.0) - args(1) * args(1))),
               static_cast<fp_type>(2.0));
  auto term_3 = std::pow(
      (static_cast<fp_type>(2.625) -
       args(0) * (static_cast<fp_type>(1.0) - args(1) * args(1) * args(1))),
      static_cast<fp_type>(2.0));
  return (term_1 + term_2 + term_3);
}

/**
 * @brief Wood's test function
 * @details initial guess = (-3.0,-1.0,-3.0,-1.0),
 * minimiser = (1.0,1.0,1.0,1.0)
 *
 * @tparam fp_type fp_type is a floation-point template parameter
 * @param args function arguments
 * @return fp_type
 */
template <typename fp_type>
fp_type wood_function(vector_arg_t<fp_type> const &args) {
  auto term_1 =
      std::pow(static_cast<fp_type>(100.0) * (args(1) - args(0) * args(0)),
               static_cast<fp_type>(2.0));
  auto term_2 = std::pow((static_cast<fp_type>(1.0) - args(0)),
                         static_cast<fp_type>(2.0));
  auto term_3 =
      std::pow(static_cast<fp_type>(90.0) * (args(3) - args(2) * args(2)),
               static_cast<fp_type>(2.0));
  auto term_4 = std::pow((static_cast<fp_type>(1.0) - args(2)),
                         static_cast<fp_type>(2.0));
  auto term_5 = std::pow(static_cast<fp_type>(10.0) *
                             (args(1) + args(3) - static_cast<fp_type>(2.0)),
                         static_cast<fp_type>(2.0));
  auto term_6 = std::pow(static_cast<fp_type>(0.1) * (args(1) - args(3)),
                         static_cast<fp_type>(2.0));
  return (term_1 + term_2 + term_3 + term_4 + term_5 + term_6);
}

/**
 * @brief Create a rao test collection object
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 * @return std::vector<sptr_t<minimizer_helper<fp_type>>>
 */
template <typename fp_type>
std::vector<sptr_t<minimizer_helper<fp_type>>> create_rao_test_collection() {

  std::vector<sptr_t<minimizer_helper<fp_type>>> helper;

  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Rosenbrock's parabolic valley"s, rosenbrock_parabolic_valley<fp_type>,
      vector_const_t<2, fp_type>(-1.2, 1.0),
      vector_const_t<2, fp_type>(1.0, 1.0)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Quadratic function"s, quadratic_function<fp_type>,
      vector_const_t<2, fp_type>(0.0, 0.0),
      vector_const_t<2, fp_type>(1.0, 3.0)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Powell's quadratic function"s, powell_function<fp_type>,
      vector_const_t<4, fp_type>(3.0, -1.0, 0.0, 1.0),
      vector_const_t<4, fp_type>(0.0, 0.0, 0.0, 0.0)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Fletcher and Powell's helical valley"s,
      fletcher_powell_helical_valley<fp_type>,
      vector_const_t<4, fp_type>(-1.0, 0.0, 0.0, 0.0),
      vector_const_t<4, fp_type>(1.0, 0.0, 0.0, 0.0)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Non-linear function of 3 variables"s, non_linear_function<fp_type>,
      vector_const_t<3, fp_type>(0.0, 1.0, 2.0),
      vector_const_t<3, fp_type>(1.0, 1.0, 1.0)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Beale's function"s, beale_function<fp_type>,
      vector_const_t<2, fp_type>(1.0, 1.0),
      vector_const_t<2, fp_type>(3.0, 0.5)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Freudenstein and Roth function"s, freudenstein_roth_function<fp_type>,
      vector_const_t<2, fp_type>(0.5, -2.0),
      vector_const_t<2, fp_type>(11.41, -0.8968)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Powell's badly scaled function"s, powell_badly_scaled_function<fp_type>,
      vector_const_t<2, fp_type>(0.0, 1.0),
      vector_const_t<2, fp_type>(1.098e-5, 9.106)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Wood's function"s, wood_function<fp_type>,
      vector_const_t<4, fp_type>(-3.0, -1.0, -3.0, -1.0),
      vector_const_t<4, fp_type>(1.0, 1.0, 1.0, 1.0)));

  return helper;
}

} // namespace om_test_functions
#endif /// OM_TEST_FUNCTIONS
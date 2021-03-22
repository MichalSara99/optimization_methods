#if !defined(OM_RAOS_COLLECTION)
#define OM_RAOS_COLLECTION

#include "om_test_helpers.hpp"
#include "utilities/om_types.hpp"
#include <cmath>
#include <memory>

using namespace om_types;
using namespace om_test_helpers;
using namespace std::string_literals;

/// Some classical test functions (designed by Rao) are listed below:
template <typename fp_type> constexpr fp_type pi{3.14159265359};

// Rosenbrock's parabolic valley
// guess = (-1.2,1.0)
// minimiser = (1.0,1.0)
template <typename fp_type>
fp_type rosenbrock_parabolic_valley(vector_arg_t<fp_type> const &args) {
  return (static_cast<fp_type>(100.0) * (args(1) - args(0) * args(0)) *
              (args(1) - args(0) * args(0)) +
          (static_cast<fp_type>(1.0) - args(0)) *
              (static_cast<fp_type>(1.0) - args(0)));
}

// Quadratic function
// guess = (0.0,0.0)
// minimiser = (1.0,3.0)
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

// Powell's quadratic function
// guess = (3.0,-1.0,0.0,1.0)
// minimiser = (0.0,0.0,0.0,0.0)
template <typename fp_type>
fp_type powell_function(vector_arg_t<fp_type> const &args) {
  return ((args(0) + static_cast<fp_type>(10.0) * args(1)) *
              (args(0) + static_cast<fp_type>(10.0) * args(1)) +
          5.0 * (args(2) - args(3)) * (args(2) - args(3)) +
          std::pow((args(1) - static_cast<fp_type>(2.0) * args(2)),
                   static_cast<fp_type>(4.0)) +
          std::pow(static_cast<fp_type>(10.0) * (args(0) - args(3)),
                   static_cast<fp_type>(4.0)));
}

// Fletcher and Powell's helical valley
// guess = (-1.0,0.0,0.0)
// minimiser = (1.0,0.0,0.0)
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

// Non-linear function of 3 variables
// guess = (0.0,1.0,2.0)
// minimiser = (1.0,1.0,1.0)
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

// Freudenstein and Roth function
// guess = (0.5,-2.0)
// minimiser = (5.0,4.0)
// local_minimiser = (11.41..., -0.8968)
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

// Powell's badly scaled function
// guess = (0.0,1.0)
// minimiser = (1.098...*10^-5,9.106...)
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

// Beale's function
// guess = (1.0,1.0)
// minimiser = (3.0,0.5)
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

// Wood's function
// guess = (-3.0,-1.0,-3.0,-1.0)
// minimiser = (1.0,1.0,1.0,1.0)
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

template <typename fp_type>
std::vector<sptr_t<minimizer_helper<fp_type>>> create_rao_test_collection() {

  std::vector<sptr_t<minimizer_helper<fp_type>>> helper;

  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Rosenbrock's parabolic valley"s, rosenbrock_parabolic_valley<fp_type>,
      Eigen::Vector2d(-1.2, 1.0), Eigen::Vector2d(1.0, 1.0)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Quadratic function"s, quadratic_function<fp_type>,
      Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(1.0, 3.0)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Powell's quadratic function"s, powell_function<fp_type>,
      Eigen::Vector4d(3.0, -1.0, 0.0, 1.0),
      Eigen::Vector4d(0.0, 0.0, 0.0, 0.0)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Fletcher and Powell's helical valley"s,
      fletcher_powell_helical_valley<fp_type>,
      Eigen::Vector4d(-1.0, 0.0, 0.0, 0.0),
      Eigen::Vector4d(1.0, 0.0, 0.0, 0.0)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Non-linear function of 3 variables"s, non_linear_function<fp_type>,
      Eigen::Vector3d(0.0, 1.0, 2.0), Eigen::Vector3d(1.0, 1.0, 1.0)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Freudenstein and Roth function"s, freudenstein_roth_function<fp_type>,
      Eigen::Vector2d(0.5, -2.0), Eigen::Vector2d(11.41, -0.8968)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Powell's badly scaled function"s, powell_badly_scaled_function<fp_type>,
      Eigen::Vector2d(0.0, 1.0), Eigen::Vector2d(1.098e-5, 9.106)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Beale's function"s, beale_function<fp_type>, Eigen::Vector2d(1.0, 1.0),
      Eigen::Vector2d(3.0, 0.5)));
  helper.emplace_back(std::make_shared<minimizer_helper<fp_type>>(
      "Wood's function"s, wood_function<fp_type>,
      Eigen::Vector4d(-3.0, -1.0, -3.0, -1.0),
      Eigen::Vector4d(1.0, 1.0, 1.0, 1.0)));

  return helper;
}

#endif /// OM_RAOS_COLLECTION
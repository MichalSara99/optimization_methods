#if !defined(OM_BRENT)
#define OM_BRENT

#include <limits>
#include <type_traits>

#include "utilities/om_types.hpp"
#include "utilities/om_utilities.hpp"

namespace om_unconstrained_methods {

namespace om_line_methods {

using om_types::f_scalar_t;
using om_utilities::fib;
using om_utilities::iqerp;
using om_utilities::lerp;
using om_utilities::range;
using om_utilities::sign;

template <typename fp_type = double,
          typename = typename std::enable_if<
              std::is_floating_point<fp_type>::value>::type>
class brent_method {
private:
  std::size_t max_iters_;
  range<fp_type> range_;
  fp_type tol_;

public:
  typedef fp_type value_type;

  brent_method(range<fp_type> const &range, fp_type tolerance = 1e-5,
               std::size_t max_iters = 1000)
      : range_{range}, tol_{tolerance}, max_iters_{max_iters} {}

  virtual ~brent_method() {}

  brent_method(brent_method const &copy)
      : range_{copy.range_}, tol_{copy.tol_}, max_iters_{copy.max_iters_} {}

  brent_method &operator=(brent_method const &copy) {
    if (&copy != this) {
      range_ = copy.range_;
      tol_ = copy.tol_;
      max_iters_ = copy.max_iters_;
    }
    return *this;
  }

  std::tuple<fp_type, fp_type, std::size_t, std::size_t>
  operator()(f_scalar_t<fp_type> &&fun) const {

    fp_type eps_2 =
        static_cast<fp_type>(2.0) * std::numeric_limits<fp_type>::epsilon();
    fp_type x0{range_.low()};
    fp_type x1{range_.high()};
    fp_type y0 = fun(x0);
    fp_type y1 = fun(x1);
    fp_type tmpx{};
    fp_type tmpy{};
    bool biFlag{true};
    // Swapping lower and upper bounds:
    if (std::abs(y0) < std::abs(y1)) {
      tmpx = x0;
      tmpy = y0;
      x0 = x1;
      x1 = tmpx;
      y0 = y1;
      y1 = tmpy;
    }

    fp_type x2{x0};
    fp_type y2{y0};
    fp_type x3{x2};
    fp_type newX{};
    fp_type newY{};
    fp_type dx{};
    fp_type dxx1{};
    fp_type dx1x2{};
    fp_type dx2x3{};
    std::size_t fun_evals{2};
    std::size_t iters{0};

    while (iters < max_iters_) {
      iters += 1;
      // if x0 and x1 are close enough we are done:
      if (std::abs(x0 - x1) < tol_)
        return std::make_tuple(x0, y0, fun_evals, iters);

      // Inverse quadratic interpolation if y0 <> y1 <> y2 else
      // secant method (linear interpolation)
      if ((std::abs(y0 - y2) > eps_2) && (std::abs(y1 - y2) > eps_2))
        newX = iqerp(x0, x1, x2, y0, y1, y2);
      else
        newX = lerp(x0, x1, y0, y1);

      // Backup bisection method if conditions satisfied:
      dx = std::abs(eps_2 * std::abs(x1));
      dxx1 = std::abs(newX - x1);
      dx1x2 = std::abs(x1 - x2);
      dx2x3 = std::abs(x2 - x3);

      if ((((static_cast<fp_type>(3.0) * x0 + x1) / static_cast<fp_type>(4.0)) <
           newX) &&
              (newX > x1) ||
          (biFlag && (dxx1 >= (static_cast<fp_type>(0.5) * dx1x2))) ||
          (!biFlag && (dxx1 >= (static_cast<fp_type>(0.5) * dx2x3))) ||
          (biFlag && (dx1x2 < dx)) || (!biFlag && (dx2x3 < dx))) {
        newX = static_cast<fp_type>(0.5) * (x0 + x1);
        biFlag = true;
      } else {
        biFlag = false;
      }

      fun_evals++;
      newY = fun(newX);

      if (std::abs(newY) < eps_2)
        return std::make_tuple(newX, newY, fun_evals, iters);

      x3 = x2;
      x2 = x1;
      if (sign(y0) != sign(newY)) {
        x1 = newX;
        y1 = newY;
      } else {
        x0 = newX;
        y0 = newY;
      }

      // Swapping lower and upper bounds:
      if (std::abs(y0) < std::abs(y1)) {
        tmpx = x0;
        tmpy = y0;
        x0 = x1;
        x1 = tmpx;
        y0 = y1;
        y1 = tmpy;
      }
    }
    return std::make_tuple(static_cast<fp_type>(0.5) * (x0 + x1),
                           static_cast<fp_type>(0.5) * (y0 + y1), fun_evals,
                           iters);
  }
};

} // namespace om_line_methods
} // namespace om_unconstrained_methods
#endif /// OM_BRENT
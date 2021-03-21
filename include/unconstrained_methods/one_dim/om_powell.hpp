#if !defined(OM_POWELL)
#define OM_POWELL

#include <limits>
#include <type_traits>

#include "utilities/om_common.hpp"
#include "utilities/om_differentiation.hpp"
#include "utilities/om_types.hpp"
#include "utilities/om_utilities.hpp"

namespace om_line_methods {

using om_common::closest_to;
using om_common::furthest_from;
using om_common::max_arg;
using om_common::min_arg;
using om_differentiation::divided_difference;
using om_types::f_scalar_t;
using om_utilities::range;

template <typename fp_type = double,
          typename = typename std::enable_if<
              std::is_floating_point<fp_type>::value>::type>
class powell_method {
private:
  range<fp_type> range_;
  fp_type tol_;
  fp_type step_;
  fp_type max_step_;
  fp_type const step_prc_{static_cast<fp_type>(0.01)};
  fp_type const max_step_prc_{static_cast<fp_type>(2.0) * step_prc_};
  std::size_t max_iters_;

public:
  typedef fp_type value_type;

  powell_method(range<fp_type> const &range, fp_type tolerance = 1e-5,
                std::size_t max_ites = 1000)
      : range_{range}, max_iters_{max_ites}, tol_{tolerance} {
    step_ = step_prc_ * range.spread();
    max_step_ = max_step_prc_ * range.spread();
  }

  powell_method(range<fp_type> const &range, fp_type step, fp_type max_step,
                fp_type tolerance = 1e-5, std::size_t max_ites = 1000)
      : range_{range}, step_{step}, max_step_{max_step},
        max_iters_{max_ites}, tol_{tolerance} {}

  virtual ~powell_method() {}

  powell_method(powell_method const &copy)
      : range_{copy.range_}, tol_{copy.tol_}, step_{copy.step_},
        max_step_{copy.max_step_} {}

  powell_method &operator=(powell_method const &copy) {
    if (this != &copy) {
      range_ = copy.range_;
      tol_ = copy.tol_;
      step_ = copy.step_;
      max_step_ = copy.max_step_;
    }
    return *this;
  }

  std::tuple<fp_type, fp_type, std::size_t, std::size_t>
  operator()(f_scalar_t<fp_type> &&fun) const {

    fp_type low = range_.low();
    fp_type high = range_.high();

    fp_type x0 = low;
    fp_type y0 = fun(x0);
    fp_type x1 = low + step_;
    fp_type y1 = fun(x1);
    fp_type x2{0.0};
    fp_type y2{0.0};

    fp_type con{0.0};

    fp_type x{0.0};
    fp_type y{0.0};

    std::size_t fun_evals{2};
    std::size_t iters{0};

    if (y0 < y1) {
      x2 = low - step_;
      y2 = fun(x2);
    } else {
      x2 = low + static_cast<fp_type>(2.0) * step_;
      y2 = fun(x2);
    }

    std::pair<fp_type, fp_type> xmin;
    std::pair<fp_type, fp_type> xmax;
    std::pair<fp_type, fp_type> xmax_2;
    std::pair<fp_type, fp_type> xoptim;
    fp_type xclosest{0.0};
    fp_type xfurthest{0.0};

    while (iters < max_iters_) {
      xmin = min_arg<3>()(fun, x0, x1, x2); // funEvals: 6
      xmax = max_arg<3>()(fun, x0, x1, x2); // funEvals: 6

      x = ((divided_difference<2>()(fun, {x0, x1, x2}) * (x0 + x1) -
            divided_difference<1>()(fun, {x0, x1})) /
           (static_cast<fp_type>(2.0) *
            divided_difference<2>()(fun, {x0, x1, x2}))); // funEvel: 4+4+2 = 10

      con = divided_difference<2>()(fun, {x0, x1, x2}); // funEval: 4
      xclosest = closest_to<3>()(x, x0, x1, x2);
      xfurthest = furthest_from<3>()(x, x0, x1, x2);

      fun_evals += 26;

      if ((con > 0) && (std::abs(x - xclosest) >
                        max_step_)) { // x is minimum | this branch: funEval: 2
                                      // discard furtest point from x:
        auto f_minus = fun(xmin.first - max_step_);
        auto f_plus = fun(xmin.first + max_step_);
        if (x0 == xfurthest) {
          if (f_minus < f_plus) {
            x0 = xmin.first - max_step_;
          } else {
            x0 = xmin.first + max_step_;
          }
        } else if (x1 == xfurthest) {
          if (f_minus < f_plus) {
            x1 = xmin.first - max_step_;
          } else {
            x1 = xmin.first + max_step_;
          }
        } else {
          if (f_minus < f_plus) {
            x2 = xmin.first - max_step_;
          } else {
            x2 = xmin.first + max_step_;
          }
        }
        fun_evals += 2;

      } else if (con <= 0) { // x is maximum | this branch funEval: 2
                             // discard closest point to x
        auto f_minus = fun(xmax.first - max_step_);
        auto f_plus = fun(xmax.first + max_step_);
        if (x0 == xclosest) {
          if (f_minus < f_plus) {
            x0 = xmax.first - max_step_;
          } else {
            x0 = xmax.first + max_step_;
          }
        } else if (x1 == xclosest) {
          if (f_minus < f_plus) {
            x1 = xmax.first - max_step_;
          } else {
            x1 = xmax.first + max_step_;
          }
        } else {
          if (f_minus < f_plus) {
            x2 = xmax.first - max_step_;
          } else {
            x2 = xmax.first + max_step_;
          }
        }
        fun_evals += 2;
      } else if (std::abs(x - xclosest) < tol_) {
        xoptim = min_arg<2>()(fun, x, xclosest); // funEval: 2
        fun_evals += 2;
        break;
      } else {
        if (x0 == xmax.first) {
          if (((x0 < x) && (x >= x1) && (x >= x2)) ||
              ((x1 <= x) && (x2 <= x) && (x < x0))) {
            xmax_2 = max_arg<2>()(fun, x1, x2);
            if (x1 == xmax_2.first) {
              x1 = x;
            } else {
              x2 = x;
            }
          } else {
            x0 = x;
          }
        } else if (x1 == xmax.first) {
          if (((x1 < x) && (x >= x0) && (x >= x2)) ||
              ((x0 <= x) && (x2 <= x) && (x < x1))) {
            xmax_2 = max_arg<2>()(fun, x0, x2);
            if (x0 == xmax_2.first) {
              x0 = x;
            } else {
              x2 = x;
            }
          } else {
            x1 = x;
          }
        } else {
          if (((x2 < x) && (x >= x0) && (x >= x1)) ||
              ((x0 <= x) && (x1 <= x) && (x < x2))) {
            xmax_2 = max_arg<2>()(fun, x0, x1);
            if (x0 == xmax_2.first) {
              x0 = x;
            } else {
              x1 = x;
            }
          } else {
            x2 = x;
          }
        }
      }
      iters++;
    }
    return std::make_tuple(xoptim.first, xoptim.second, fun_evals, iters);
  }
};

} // namespace om_line_methods

#endif // OM_POWELL
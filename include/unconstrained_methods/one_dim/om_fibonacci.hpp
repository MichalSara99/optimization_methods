#if !defined(OM_FIBONACCI)
#define OM_FIBONACCI

#include "utilities/om_types.hpp"
#include "utilities/om_utilities.hpp"

#include <type_traits>

namespace om_unconstrained_methods {

namespace om_line_methods {

using om_types::f_scalar_t;
using om_utilities::fib;
using om_utilities::range;

/**
 * @brief Fibonacci method object
 *
 * @tparam fp_type fp_type is floating-point template parameter
 * @tparam std::enable_if<
 * std::is_floating_point<fp_type>::value>::type
 */
template <typename fp_type = double,
          typename = typename std::enable_if<
              std::is_floating_point<fp_type>::value>::type>
class fibonacci_method {
private:
  std::size_t max_iters_;
  range<fp_type> range_;
  fp_type tol_;

public:
  typedef fp_type value_type;
  /**
   * @brief Construct a new fibonacci method object
   *
   *
   * @param range range of the minimiser
   * @param tolerance tolerance of the minimiser
   * @param max_iters maximum number of iterations
   */
  fibonacci_method(range<fp_type> const &range, fp_type tolerance = 1e-5,
                   std::size_t max_iters = 1000)
      : range_{range}, tol_{tolerance}, max_iters_{max_iters} {}

  virtual ~fibonacci_method() {}
  /**
   * @brief Copy constructor of a fibonacci method object
   *
   * @param copy copy is the object which we want to make a copy of
   */
  fibonacci_method(fibonacci_method const &copy)
      : range_{copy.range_}, tol_{copy.tol_}, max_iters_{copy.max_iters_} {}
  /**
   * @brief Assignment operator of a fibonacci method object
   *
   * @param copy
   * @return fibonacci_method&
   */
  fibonacci_method &operator=(fibonacci_method const &copy) {
    if (&copy != this) {
      range_ = copy.range_;
      tol_ = copy.tol_;
      max_iters_ = copy.max_iters_;
    }
    return *this;
  }
  /**
   * @brief Functor of a fibonacci method object
   *
   *
   * @param fun objective function
   * @return std::tuple<fp_type, fp_type, std::size_t, std::size_t>
   */
  std::tuple<fp_type, fp_type, std::size_t, std::size_t>
  operator()(f_scalar_t<fp_type> &&fun) const {

    std::size_t N{max_iters_};
    fp_type r{static_cast<fp_type>(fib(N) / fib(N + 1))};
    fp_type L{range_.spread()};
    fp_type low{range_.low()};
    fp_type high{range_.high()};
    fp_type lambda_1{low + r * r * L};
    fp_type lambda_2{low + r * L};
    std::size_t fun_evals{0};
    std::size_t iters{0};

    while ((L >= tol_) && (iters < max_iters_)) {
      iters += 1;
      r = (fib(N - iters + 1) / fib(N - iters + 2));
      if (fun(lambda_1) > fun(lambda_2)) {
        low = lambda_1;
        lambda_1 = lambda_2;
        L = high - low;
        lambda_2 = low + r * L;
      } else {
        high = lambda_2;
        lambda_2 = lambda_1;
        L = high - low;
        lambda_1 = low + r * r * L;
      }
      fun_evals += 2;
    }
    return std::make_tuple(static_cast<fp_type>(0.5) * (low + high),
                           fun(static_cast<fp_type>(0.5) * (low + high)),
                           fun_evals, iters);
  }
};

} // namespace om_line_methods
} // namespace om_unconstrained_methods
#endif /// OM_FIBONACCI
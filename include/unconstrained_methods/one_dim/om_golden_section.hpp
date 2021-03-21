#if !defined(OM_GOLDEN_SECTION)
#define OM_GOLDEN_SECTION

#include "utilities/om_types.hpp"
#include "utilities/om_utilities.hpp"

#include <limits>
#include <type_traits>

namespace om_line_methods {

using om_types::f_scalar_t;
using om_utilities::range;

template <typename fp_type = double,
          typename = typename std::enable_if<
              std::is_floating_point<fp_type>::value>::type>
class golden_section_method {
private:
  std::size_t max_iters_;
  range<fp_type> range_;
  fp_type tol_;
  const fp_type golden_section_{
      static_cast<fp_type>(0.5 * (std::sqrt(5.0) - 1.0))};

public:
  typedef fp_type value_type;

  golden_section_method(range<fp_type> const &range, fp_type tolerance = 1e-5,
                        std::size_t max_iters = 1000)
      : range_{range}, tol_{tolerance}, max_iters_{max_iters} {}

  virtual ~golden_section_method() {}

  golden_section_method(golden_section_method const &copy)
      : range_{copy.range_}, tol_{copy.tol_}, max_iters_{copy.max_iters_},
        golden_section_{copy.golden_section_} {}

  golden_section_method &operator=(golden_section_method const &copy) {
    if (&copy != this) {
      range_ = copy.range_;
      tol_ = copy.tol_;
      max_iters_ = copy.max_iters_;
    }
    return *this;
  }

  std::tuple<fp_type, fp_type, std::size_t, std::size_t>
  operator()(f_scalar_t<fp_type> &&fun) const {

    fp_type const r = golden_section_;
    fp_type L{range_.spread()};
    fp_type low{range_.low()};
    fp_type high{range_.high()};
    fp_type lambda_1{low + r * r * L};
    fp_type lambda_2{low + r * L};
    std::size_t fun_evals{0};
    std::size_t iters{0};

    while ((L >= tol_) && (iters < max_iters_)) {
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
      iters += 1;
    }
    return std::make_tuple(static_cast<fp_type>(0.5) * (low + high),
                           fun(static_cast<fp_type>(0.5) * (low + high)),
                           fun_evals, iters);
  }
};
} // namespace om_line_methods

#endif /// OM_GOLDEN_SECTION
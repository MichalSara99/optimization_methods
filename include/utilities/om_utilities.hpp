#ifndef OM_UTILITIES
#define OM_UTILITIES

#include "om_types.hpp"
#include <random>
#include <vector>

/**
 * @brief Contains some commonly used utilities
 *
 */
namespace om_utilities {

using om_types::vector_t;

/**
 * @brief Represents a one dimensional range
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 * @tparam std::enable_if<
 * std::is_floating_point<fp_type>::value>::type
 */
template <typename fp_type = double,
          typename = typename std::enable_if<
              std::is_floating_point<fp_type>::value>::type>
class range {
private:
  fp_type low_;
  fp_type high_;

public:
  /**
   * @brief Construct a new range object
   *
   * @param low low value of a range
   * @param high high value of a range
   */
  range(fp_type const &low, fp_type const &high) : low_{low}, high_{high} {}
  /**
   * @brief Construct a new range object
   *
   */
  range() : range(fp_type{}, fp_type{}) {}
  /**
   * @brief Construct a new range object
   *
   * @param copy copy is the object which we want to make a copy of
   */
  range(range<fp_type> const &copy) : low_{copy.low_}, high_{copy.high_} {}

  /**
   * @brief Copy assignment operator of a range object
   *
   * @param copy
   * @return range<fp_type>&
   */
  range<fp_type> &operator=(range<fp_type> const &copy) {
    if (this != &copy) {
      low_ = copy.low_;
      high_ = copy.high_;
    }
  }
  /**
   * @brief Move constructor of a range object
   *
   * @param other
   */
  range(range<fp_type> &&other) {
    low_ = std::exchange(other.low_, fp_type{});
    high_ = std::exchange(other.high_, fp_type{});
  }

  /**
   * @brief Move assignment of a range object
   *
   * @param other
   * @return range<fp_type>&
   */
  range<fp_type> &operator=(range<fp_type> &&other) {
    if (this != &other) {
      low_ = std::exchange(other.low_, fp_type{});
      high_ = std::exchange(other.high_, fp_type{});
    }
  }
  /**
   * @brief Returns low end of the range
   *
   * @return fp_type const&
   */
  inline fp_type const &low() const { return low_; }
  /**
   * @brief Returns high end of the range
   *
   * @return fp_type const&
   */
  inline fp_type const &high() const { return high_; }
  /**
   * @brief Returns a pair of low high end of the range
   *
   * @return std::pair<fp_type, fp_type>
   */
  inline std::pair<fp_type, fp_type> low_high() const {
    return std::make_pair(low_, high_);
  }
  /**
   * @brief Returns a spread between high and low end of the range
   *
   * @return fp_type
   */
  inline fp_type spread() const { return (high_ - low_); }
};

/**
 * @brief Random vectors from guess functor
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 * @tparam distribution distribution of random generator
 * @tparam std::enable_if_t<std::is_floating_point<fp_type>::value>
 */
template <typename fp_type = double,
          template <typename> typename distribution = std::normal_distribution,
          typename =
              typename std::enable_if_t<std::is_floating_point<fp_type>::value>>
struct random_vectors_from_guess {
private:
  std::random_device rd_;
  distribution<fp_type> dist_;

public:
  std::vector<vector_t<fp_type>>
  operator()(std::size_t N, vector_t<fp_type> const &init_guess) {
    auto mtgen = std::mt19937{rd_()};
    std::vector<vector_t<fp_type>> container(N);

    for (std::size_t t = 0; t < N; ++t) {
      auto new_vector = init_guess;
      for (auto r = 0; r < init_guess.rows(); ++r) {
        new_vector(r) += dist_(mtgen);
      }
      container[t] = std::move(new_vector);
    }
    return container;
  }
};

/**
 * @brief Cartesian basis vectors functor
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 * @tparam std::enable_if<
 * std::is_floating_point<fp_type>::value>::type
 */
template <typename fp_type = double,
          typename = typename std::enable_if<
              std::is_floating_point<fp_type>::value>::type>
struct cartesian_basis_vectors {
public:
  std::vector<vector_t<fp_type>>
  operator()(std::size_t const &dimension) const {
    std::vector<vector_t<fp_type>> basis(dimension);
    for (auto t = 0; t < dimension; ++t) {
      vector_t<fp_type> e = vector_t<fp_type>::Zero(dimension);
      e(t) = static_cast<fp_type>(1.0);
      basis[t] = std::move(e);
    }
    return basis;
  }
};

/**
 * @brief fib function
 *
 * @param n number of values from Fibonacci sequence
 * @return double
 */
double fib(std::size_t n) {
  if (n == 0)
    return 1;
  if (n == 1)
    return 1;
  double f_0 = 1;
  double f_1 = 1;
  double f_n{0};
  std::size_t i{2};

  while (i <= n) {
    f_n = f_1 + f_0;
    f_0 = f_1;
    f_1 = f_n;
    i++;
  }
  return f_n;
}

/**
 * @brief Inverse quadratic interpolation among points
 * (x0,y0),(x1,y1),(x2,y2)
 *
 * @tparam fp_type
 * @param x0 first value
 * @param x1 second value
 * @param x2 third value
 * @param y0 first function value
 * @param y1 second function value
 * @param y2 third function value
 * @return fp_type
 */
template <typename fp_type>
fp_type iqerp(fp_type x0, fp_type x1, fp_type x2, fp_type y0, fp_type y1,
              fp_type y2) {
  fp_type t1 = (x0 * y1 * y2 / ((y0 - y1) * (y0 - y2)));
  fp_type t2 = (x1 * y0 * y2 / ((y1 - y0) * (y1 - y2)));
  fp_type t3 = (x2 * y0 * y1 / ((y2 - y0) * (y2 - y1)));
  return (t1 + t2 + t3);
}

/**
 * @brief Linear interpolation between points (x0,y0) and (x1,y1)
 *
 * @tparam fp_type
 * @param x0 first value
 * @param x1 second value
 * @param y0 first function value
 * @param y1 second function value
 * @return fp_type
 */
template <typename fp_type>
fp_type lerp(fp_type x0, fp_type x1, fp_type y0, fp_type y1) {
  return (x1 - y1 * (x1 - x0) / (y1 - y0));
}

/**
 * @brief Signum function
 *
 * @tparam fp_type fp_type is a floating-point template parameter
 * @param x value
 * @return fp_type
 */
template <typename fp_type> fp_type sign(fp_type x) {
  if (x < static_cast<fp_type>(0.0))
    return static_cast<fp_type>(-1.0);
  else if (x > static_cast<fp_type>(0.0))
    return static_cast<fp_type>(1.0);
  return static_cast<fp_type>(0.0);
}

} // namespace om_utilities

#endif /// OM_UTILITIES
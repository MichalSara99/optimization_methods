#ifndef OM_UTILITIES
#define OM_UTILITIES

#include "om_types.hpp"
#include <random>
#include <vector>

namespace om_utilities {

using om_types::vector_t;

// =============================================================
// ========================== range ============================
// =============================================================
template <typename fp_type = double,
          typename = typename std::enable_if<
              std::is_floating_point<fp_type>::value>::type>
class range {
private:
  fp_type low_;
  fp_type high_;

public:
  range(fp_type const &low, fp_type const &high) : low_{low}, high_{high} {}
  range() : range(fp_type{}, fp_type{}) {}

  range(range<fp_type> const &copy) : low_{copy.low_}, high_{copy.high_} {}

  range<fp_type> &operator=(range<fp_type> const &copy) {
    if (this != &copy) {
      low_ = copy.low_;
      high_ = copy.high_;
    }
  }

  range(range<fp_type> &&other) {
    low_ = std::exchange(other.low_, fp_type{});
    high_ = std::exchange(other.high_, fp_type{});
  }

  range<fp_type> &operator=(range<fp_type> &&other) {
    if (this != &other) {
      low_ = std::exchange(other.low_, fp_type{});
      high_ = std::exchange(other.high_, fp_type{});
    }
  }

  inline fp_type const &low() const { return low_; }
  inline fp_type const &high() const { return high_; }
  inline std::pair<fp_type, fp_type> low_high() const {
    return std::make_pair(low_, high_);
  }
  inline fp_type spread() const { return (high_ - low_); }
};

// ==================================================================
// =================== random_vectors_from_guess ====================
// ==================================================================

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

// ============================================================================
// ========================== cartesian_basis_vectors =========================
// ============================================================================

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

// ================================================================
// ============================== fib =============================
// ================================================================

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

// ================================================================
// ================ inverse quadratic interpolation ===============
// ================================================================

template <typename T> T iqerp(T x0, T x1, T x2, T y0, T y1, T y2) {
  T t1 = (x0 * y1 * y2 / ((y0 - y1) * (y0 - y2)));
  T t2 = (x1 * y0 * y2 / ((y1 - y0) * (y1 - y2)));
  T t3 = (x2 * y0 * y1 / ((y2 - y0) * (y2 - y1)));
  return (t1 + t2 + t3);
}

// ================================================================
// ======================= linear interpolation ===================
// ================================================================

template <typename T> T lerp(T x0, T x1, T y0, T y1) {
  return (x1 - y1 * (x1 - x0) / (y1 - y0));
}

// ================================================================
// ============================== sign ============================
// ================================================================

template <typename T> T sign(T x) {
  if (x < 0.0)
    return static_cast<T>(-1.0);
  else if (x > static_cast<T>(0.0))
    return static_cast<T>(1.0);
  return static_cast<T>(0.0);
}

} // namespace om_utilities

#endif /// OM_UTILITIES
#if !defined(OM_TYPES)
#define OM_TYPES

#if __linux__
#include <Eigen/Core>
#elif _WIN32
#include <Eigen/Core.hpp>
#endif
#include <complex>
#include <functional>
#include <map>
#include <memory>
#include <vector>

namespace om_types {
template <typename T> using sptr_t = std::shared_ptr<T>;

template <typename T> using vector_arg_t = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <typename T> using vector_t = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <typename T> using f_scalar_t = std::function<T(T)>;

template <typename T> using f_vector_t = std::function<T(vector_arg_t<T>)>;

template <typename T>
using matrix_t = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename fp_type>
using f_line_minimiser_t =
    std::function<std::tuple<fp_type, fp_type, std::size_t, std::size_t>(
        f_scalar_t<fp_type> &&)>;

enum class one_dim_line_search_method { GoldenSection, Powell };

enum class constraint_t { Equality, LessThenZero };

template <typename T>
using constraints_t = std::vector<std::pair<f_vector_t<T>, constraint_t>>;
} // namespace om_types

#endif /// OM_TYPES
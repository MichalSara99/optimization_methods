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

/**
 * @brief Contains some types used throughout the whole library
 *
 */
namespace om_types {
/**
 * @brief Alias for shared_ptr<T>
 *
 * @tparam T
 */
template <typename T> using sptr_t = std::shared_ptr<T>;
/**
 * @brief Alias for 1D matrix = vector
 *
 * @tparam T
 */
template <typename T> using vector_arg_t = Eigen::Matrix<T, Eigen::Dynamic, 1>;
/**
 * @brief Alias for const dimension 1D matrix = vector<dimension>
 *
 * @tparam dimension
 * @tparam T
 */
template <std::size_t dimension, typename T>
using vector_const_t = Eigen::Matrix<T, dimension, 1>;
/**
 * @brief Alias for dynamic 1D matrix = vector
 *
 * @tparam T
 */
template <typename T> using vector_t = Eigen::Matrix<T, Eigen::Dynamic, 1>;
/**
 * @brief One dimensional scalar function
 *
 * @tparam T
 */
template <typename T> using f_scalar_t = std::function<T(T)>;
/**
 * @brief One dimensional vector function
 *
 * @tparam T
 */
template <typename T> using f_vector_t = std::function<T(vector_arg_t<T>)>;
/**
 * @brief Alias for Eigen matrix
 *
 * @tparam T
 */
template <typename T>
using matrix_t = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
/**
 * @brief Line method functor type
 *
 * @tparam fp_type
 */
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
#pragma once
#include <omp.h>
#include <vector>
#include <cmath>
#include "types.h"

/**
 * @namespace math_utils
 * @brief Auxiliary mathematical functions: dot product, norm, constants.
 */
namespace math_utils {
/**
* @brief Pi constant.
*/
constexpr long double PI = 3.14159265358979323846;

/**
* @brief Dot product of two vectors.
* @param a First vector.
* @param b Second vector.
* @return Sum of a[i]*b[i] over all elements.
* @note Uses OpenMP parallel reduction.
*/
inline DataType dot(const std::vector<DataType>& a, const std::vector<DataType>& b) {
    DataType sum = 0.0;
#pragma omp parallel for reduction(+:sum)
    for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(a.size()); ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

/**
* @brief Euclidean norm of a vector.
* @param v Input vector.
* @return sqrt(dot(v,v)).
*/
inline DataType norm(const std::vector<DataType>& v) {
    return std::sqrt(dot(v, v));
}

}  // namespace math_utils

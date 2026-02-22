#pragma once
#include <cmath>
#include <tuple>
#include "types.h"
#include "math_utils.h"

/**
 * @namespace eigen_theory
 * @brief Theoretical formulas for eigenvalues of the Laplace operator.
 */
namespace eigen_theory {
/**
* @brief Computes the theoretical minimum eigenvalue.
* @param coeff_x Coefficient a.
* @param coeff_y Coefficient b.
* @param h Grid step.
* @return Minimum eigenvalue: 4(a+b)/h^2 * sin^2(pi*h/2).
*/
inline DataType lambda_min(DataType coeff_x, DataType coeff_y, DataType h) {
    DataType coeff_sum = coeff_x + coeff_y;
    DataType factor = 4.0 * coeff_sum / (h * h);
    DataType sin_term = std::sin(math_utils::PI * h / 2.0);
    return factor * sin_term * sin_term;
}

/**
* @brief Computes the theoretical maximum eigenvalue.
* @param coeff_x Coefficient a.
* @param coeff_y Coefficient b.
* @param h Grid step.
* @return Maximum eigenvalue: 4(a+b)/h^2 * cos^2(pi*h/2).
*/
inline DataType lambda_max(DataType coeff_x, DataType coeff_y, DataType h) {
    DataType coeff_sum = coeff_x + coeff_y;
    DataType factor = 4.0 * coeff_sum / (h * h);
    DataType cos_term = std::cos(math_utils::PI * h / 2.0);
    return factor * cos_term * cos_term;
}

/**
* @brief Returns both eigenvalues as a tuple.
* @param coeff_x Coefficient a.
* @param coeff_y Coefficient b.
* @param h Grid step.
* @return Tuple (lambda_max, lambda_min).
*/
inline std::tuple<DataType, DataType> lambda_max_min(DataType coeff_x, DataType coeff_y, DataType h) {
    return {
        lambda_max(coeff_x, coeff_y, h),
        lambda_min(coeff_x, coeff_y, h)
    };
}
} // namespace eigen_theory

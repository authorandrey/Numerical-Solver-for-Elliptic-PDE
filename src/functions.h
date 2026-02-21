#pragma once
#include "context.h"
#include <cmath>

/**
 * @namespace exact
 * @brief Exact solution and right-hand side for the test equation.
 */
namespace exact {
    /**
     * @brief Exact solution u(x,y) = x^2 + cos^2(y*x).
     * @param x x-coordinate.
     * @param y y-coordinate.
     * @return u(x,y).
     */
    inline DataType u(DataType x, DataType y) {
        return x * x + std::cos(y * x) * std::cos(y * x);
    }

    /**
     * @brief Right-hand side corresponding to the exact solution.
     * @param x x-coordinate.
     * @param y y-coordinate.
     * @param coeff_x Coefficient a.
     * @param coeff_y Coefficient b.
     * @return f(x,y) = -2a + 2 cos(2yx)(a y^2 + b x^2).
     */
    inline DataType source(DataType x, DataType y, DataType coeff_x, DataType coeff_y) {
        DataType yx = y * x;
        DataType cos2yx = std::cos(2.0 * yx);
        return -2.0 * coeff_x + 2.0 * cos2yx * (coeff_x * y * y + coeff_y * x * x);
    }

    /**
     * @brief Boundary condition (same as exact solution).
     * @param x x-coordinate.
     * @param y y-coordinate.
     * @return u(x,y).
     */
    inline DataType boundary(DataType x, DataType y) {
        return u(x, y);
    }
}
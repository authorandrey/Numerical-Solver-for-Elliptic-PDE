#pragma once
#include "mesh.h"

/**
 * @struct Context
 * @brief Contains all solver configuration parameters.
 */
struct Context {
    UniformMesh mesh;               ///< The grid.
    const DataType coeff_x;         ///< Coefficient a for -u_xx.
    const DataType coeff_y;         ///< Coefficient b for -u_yy.
    const DataType h;               ///< Grid step (duplicated from mesh).
    const DataType tol_power;       ///< Tolerance for the power method.
    const DataType tol_seidel;      ///< Tolerance for the Seidel method.
    const size_t max_iter_power;    ///< Maximum number of iterations for the power method.
    const size_t max_iter_seidel;   ///< Maximum number of iterations for the Seidel method.
    const std::string& output_dir;   ///< Directory for output files (solution, plots).
};

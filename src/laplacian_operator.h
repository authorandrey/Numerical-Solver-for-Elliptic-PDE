#pragma once
#include <omp.h>
#include <vector>
#include "mesh.h"

/**
 * @class LaplacianOperator
 * @brief Discrete operator: -a u_xx - b u_yy.
 * @details Applies the operator to a vector representing function values on the grid.
 */
class LaplacianOperator {
public:
    /**
     * @brief Constructor.
     * @param mesh Grid.
     * @param coeff_x Coefficient a.
     * @param coeff_y Coefficient b.
     */
    LaplacianOperator(const UniformMesh& mesh, DataType coeff_x, DataType coeff_y)
        : m_mesh(mesh), m_coeff_x(coeff_x), m_coeff_y(coeff_y) {
    }

    /**
     * @brief Applies the operator to vector x, result in y.
     * @param[in] x Input vector of size total_nodes().
     * @param[out] y Output vector of the same size.
     * @note Only interior nodes are computed; boundary values in x are ignored
     *       (for eigenvalue problems they should be zero).
     */
    void apply(const std::vector<DataType>& x, std::vector<DataType>& y) const {
        const DataType inv_h2 = 1.0 / (m_mesh.step() * m_mesh.step());
        const std::size_t n = m_mesh.nodes_per_side();
        const std::size_t stride = n;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (std::ptrdiff_t j = 1; j < static_cast<std::ptrdiff_t>(n - 1); ++j) {
            for (std::size_t i = 1; i < n - 1; ++i) {
                std::size_t idx = i + j * stride;
                DataType center = x[idx];
                DataType left = x[idx - 1];
                DataType right = x[idx + 1];
                DataType down = x[idx - stride];
                DataType up = x[idx + stride];

                DataType d2x = -m_coeff_x * inv_h2 * (left - 2.0 * center + right);
                DataType d2y = -m_coeff_y * inv_h2 * (down - 2.0 * center + up);

                y[idx] = d2x + d2y;
            }
        }
    }

private:
    const UniformMesh& m_mesh;
    DataType m_coeff_x, m_coeff_y;
};

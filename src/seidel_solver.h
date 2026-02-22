#pragma once
#include <omp.h>
#include <vector>
#include <cmath>
#include "mesh.h"

/**
 * @class SeidelSolver
 * @brief Solves the equation -a u_xx - b u_yy = f using the Seidel method.
 */
class SeidelSolver {
public:
    /**
     * @brief Constructor, initializes the right-hand side and boundary conditions.
     * @param mesh Grid.
     * @param coeff_x Coefficient a.
     * @param coeff_y Coefficient b.
     */
    SeidelSolver(const UniformMesh& mesh, DataType coeff_x, DataType coeff_y)
        : m_mesh(mesh), m_coeff_x(coeff_x), m_coeff_y(coeff_y),
        m_right_hand_side(m_mesh.total_nodes(), 0.0),
        m_solution(m_mesh.total_nodes(), 0.0) {
        init();
    }

    /**
     * @brief Performs Seidel iterations.
     * @param tolerance Tolerance for residual reduction.
     * @param max_iter Maximum number of iterations.
     * @param check_period Period for residual check (default 50).
     * @return Solution vector (size total_nodes()).
     * @details The residual is computed as sqrt(sum (Lu - f)^2) over interior nodes.
     */
    std::vector<DataType> solve(DataType tolerance, size_t max_iter, size_t check_period = 50) {
        const DataType h2 = m_mesh.step() * m_mesh.step();
        const DataType coeff_sum = 2.0 * (m_coeff_x + m_coeff_y);
        const size_t stride = m_mesh.nodes_per_side();

        auto residual_norm = [&]() -> DataType {
            DataType res = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:res)
#endif
            for (std::ptrdiff_t j = 1; j < static_cast<std::ptrdiff_t>(stride - 1); ++j) {
                for (size_t i = 1; i < stride - 1; ++i) {
                    size_t idx = i + j * stride;
                    DataType center = m_solution[idx];
                    DataType left = m_solution[idx - 1];
                    DataType right = m_solution[idx + 1];
                    DataType down = m_solution[idx - stride];
                    DataType up = m_solution[idx + stride];

                    DataType d2x = -m_coeff_x * (left - 2.0 * center + right) / h2;
                    DataType d2y = -m_coeff_y * (down - 2.0 * center + up) / h2;
                    DataType Lu = d2x + d2y;
                    DataType r = Lu - m_right_hand_side[idx];
                    res += r * r;
                }
            }
            return std::sqrt(res);
        };

        DataType initial_residual = residual_norm();
        DataType target_residual = tolerance * initial_residual;

        for (size_t iter = 0; iter < max_iter; ++iter) {
            for (size_t j = 1; j < stride - 1; ++j) {
                for (size_t i = 1; i < stride - 1; ++i) {
                    size_t idx = i + j * stride;
                    DataType numerator =
                        m_coeff_x * (m_solution[idx - 1] + m_solution[idx + 1]) +
                        m_coeff_y * (m_solution[idx - stride] + m_solution[idx + stride]) +
                        m_right_hand_side[idx] * h2;
                    m_solution[idx] = numerator / coeff_sum;
                }
            }

            if (iter % check_period == 0) {
                DataType current_residual = residual_norm();
                if (current_residual < target_residual)
                    break;
            }
        }

        return m_solution;
    }

private:
    /**
     * @brief Initializes the right-hand side from the exact source and sets boundary conditions.
     */
    void init() {
        for (size_t j = 0; j < m_mesh.nodes_per_side(); ++j) {
            for (size_t i = 0; i < m_mesh.nodes_per_side(); ++i) {
                DataType x = m_mesh.x(i);
                DataType y = m_mesh.y(j);
                m_right_hand_side[m_mesh.full_index(i, j)]
                    = exact::source(x, y, m_coeff_x, m_coeff_y);
                if (!m_mesh.is_inner(i, j)) {
                    m_solution[m_mesh.full_index(i, j)] = exact::boundary(x, y);
                }
            }
        }
    }

private:
    const UniformMesh& m_mesh;
    DataType m_coeff_x, m_coeff_y;
    std::vector<DataType> m_right_hand_side;
    std::vector<DataType> m_solution;
};

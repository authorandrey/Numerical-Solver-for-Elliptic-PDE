#pragma once
#include <vector>
#include <tuple>
#include "math_utils.h"
#include "context.h"

/**
 * @class PowerMethod
 * @brief Implements the power method for computing the largest eigenvalue.
 */
class PowerMethod {
public:
    /**
     * @struct Params
     * @brief Parameters for the power method.
     * @tparam Operator Type of linear operator (must have apply(x,y) method).
     */
    template<typename Operator>
    struct Params {
        Operator& op;                 ///< Reference to the operator.
        const UniformMesh& mesh;      ///< Grid for handling boundaries.
        std::vector<DataType>& x;     ///< Initial vector (modified during the process).
        const DataType tolerance;     ///< Convergence tolerance.
        const size_t max_iter;        ///< Maximum number of iterations.
    };

    /**
     * @brief Computes the largest eigenvalue using the power method.
     * @param params Parameters, including the initial vector and operator.
     * @return Tuple (eigenvalue, number of iterations).
     * @details The initial vector is normalized, then iterations:
     *          x_{k+1} = A x_k / ||A x_k||, lambda = x_k^T A x_k.
     *          Stopping criterion: residual ||A x - lambda x|| / |lambda| < epsilon.
     */
    template<typename Operator>
    std::tuple<DataType, unsigned> compute_max(Params<Operator>& params) const {
        const size_t N = params.x.size();

        normalize_with_boundary(params.x, params.mesh);

        DataType lambda = 0.0;
        std::vector<DataType> y(N, 0.0);
        unsigned iterations = 0;

        for (size_t iter = 0; iter < params.max_iter; ++iter) {
            ++iterations;
            params.op.apply(params.x, y);

            lambda = math_utils::dot(params.x, y);

            DataType residual = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:residual)
#endif
            for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(N); ++i) {
                DataType r = y[i] - lambda * params.x[i];
                residual += r * r;
            }
            residual = std::sqrt(residual);

            if (std::abs(lambda) > 0
                && residual / std::abs(lambda) < params.tolerance) {
                break;
            }

            DataType y_norm = math_utils::norm(y);
            if (y_norm == 0) break;
            // x_{k+1} = x_k / ||x_k||
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(N); ++i) {
                params.x[i] = y[i] / y_norm;
            }
            enforce_zero_boundary(params.x, params.mesh);
        }
        return { lambda, iterations };
    }

private:
    /**
     * @brief Normalizes the vector on interior nodes and sets boundaries to zero.
     * @param v Vector.
     * @param mesh Grid.
     * @details If the norm is zero, all interior nodes are filled with 1/sqrt(inner_nodes).
     */
    static void normalize_with_boundary(std::vector<DataType>& v, const UniformMesh& mesh) {
        DataType n = 0.0;
        mesh.for_each_inner(
            [&](std::size_t, std::size_t, std::size_t idx) {
                n += v[idx] * v[idx];
            }
        );
        n = std::sqrt(n);
        if (n > 0) {
            mesh.for_each_inner(
                [&](std::size_t, std::size_t, std::size_t idx) {
                    v[idx] /= n;
                }
            );
        }
        else {
            DataType val = 1.0 / std::sqrt(mesh.inner_nodes());
            mesh.for_each_inner(
                [&](std::size_t, std::size_t, std::size_t idx) {
                    v[idx] = val;
                }
            );
        }
        enforce_zero_boundary(v, mesh);
    }

    /**
     * @brief Sets all boundary nodes of the vector to zero.
     * @param v Vector.
     * @param mesh Grid.
     */
    static void enforce_zero_boundary(std::vector<DataType>& v, const UniformMesh& mesh) {
        std::size_t n = mesh.nodes_per_side();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(n); ++j) {
            v[mesh.full_index(0, j)] = 0.0;
            v[mesh.full_index(n - 1, j)] = 0.0;
        }
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (std::ptrdiff_t i = 1; i < static_cast<std::ptrdiff_t>(n - 1); ++i) {
            v[mesh.full_index(i, 0)] = 0.0;
            v[mesh.full_index(i, n - 1)] = 0.0;
        }
    }
};

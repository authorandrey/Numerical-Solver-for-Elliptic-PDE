#pragma once
#include "solvers/laplacian_operator.h"
#include "solvers/power_method.h"

/**
 * @class ShiftedPowerMethod
 * @brief Computes the smallest eigenvalue using the shifted power method.
 * @details Applies the power method to the operator shift*I - A.
 */
class ShiftedPowerMethod {
public:
    /**
     * @brief Synonim for Params in PowerMethod.
     */
    template<typename Operator>
    using Params = PowerMethod::Params<Operator>;

    /**
    * @brief Constructor with a given shift (usually the largest eigenvalue).
    * @param _shift Shift value.
    */
    explicit ShiftedPowerMethod(DataType _shift) : m_shift(_shift) { }

    /**
     * @brief Computes the smallest eigenvalue.
     * @param params Parameters for the power method (initial vector, etc.).
     * @return Tuple (eigenvalue, number of iterations).
     * @details Uses the shifted operator (shift*I - A) and then
     *          lambda_min = shift - mu, where mu is the largest eigenvalue of the shifted operator.
     */
    template<typename Operator>
    std::tuple<DataType, unsigned> compute_min(Params<Operator>& params) const {
        PowerMethod power_method;

        /**
         * @brief Internal helper class implementing the shifted operator.
         */
        class ShiftedOperator {
        public:
            ShiftedOperator(const LaplacianOperator& op, DataType shift)
                : m_op(op), m_shift(shift) {
            }

            /**
             * @brief Applies the shifted operator.
             */
            void apply(const std::vector<DataType>& x, std::vector<DataType>& y) const {
                m_op.apply(x, y);
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(x.size()); ++i) {
                    y[i] = m_shift * x[i] - y[i];
                }
            }

        private:
            const LaplacianOperator& m_op;
            DataType m_shift;
        };

        ShiftedOperator shifted_op(params.op, m_shift);

        Params shifted_params(shifted_op, params.mesh, params.x, params.tolerance, params.max_iter);

        auto [mu, iterations] = power_method.compute_max(shifted_params);
        return { m_shift - mu, iterations };
    }

private:
    const DataType m_shift;
};

#include "controller/solver_controller.h"
#include "core/settings.h"

/**
 * @brief Entry point.
 * @return 0 on success, non-zero otherwise.
 * @details Task variant: 12.
            Uses parameters: a=0.9, b=1.1, h=0.01, tolerance=1e-6, max iterations=100000.
 *          Initial vectors: Alternating for maximum, UnitConstant for minimum.
 */
int main() {
    const DataType a = 0.9, b = 1.1;
    const DataType h = 0.01;
    const DataType tol_power = 1e-6, tol_seidel = 1e-6;
    const size_t max_iter_power = 100000, max_iter_seidel = 100000;

    try {
        SolverController solver(a, b, h, tol_power, tol_seidel,
            max_iter_power, max_iter_seidel, OUTPUT_DIR);
        solver.run_all(
            InitialVectorType::Alternating,
            InitialVectorType::UnitConstant
        );
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

#pragma once
#include <random>
#include "functions.h"
#include "seidel_solver.h"
#include "power_method.h"
#include "shifted_power_method.h"
#include "initial_vector_type.h"
#include "clock_wrapper.h"
#include "eigen_theory.h"
#include "output_handler.h"
#include "context.h"

/**
 * @class SolverController
 * @brief Controller managing the solution of the equation, eigenvalue computation, and output/saving of results.
 */
class SolverController {
public:
    /**
     * @brief Constructor.
     * @param a Coefficient a.
     * @param b Coefficient b.
     * @param h Grid step.
     * @param tol_power Tolerance for the power method.
     * @param tol_seidel Tolerance for the Seidel method.
     * @param max_iter_power Maximum number of iterations for the power method.
     * @param max_iter_seidel Maximum number of iterations for the Seidel method.
     * @param output_dir Directory where output files (solution binary, plots) will be stored.
     */
    SolverController(double a, double b, double h, double tol_power, double tol_seidel,
        size_t max_iter_power, size_t max_iter_seidel, const std::string& output_dir)
        : ctx{ UniformMesh(h), a, b, h, tol_power, tol_seidel,
           max_iter_power, max_iter_seidel, output_dir } {
    }

    /**
     * @brief Runs the Seidel method, saves/plots the solution, then computes eigenvalues.
     * @param type_max Type of initial vector for the maximum eigenvalue.
     * @param type_min Type of initial vector for the minimum eigenvalue.
     */
    [[maybe_unused]]
    void run_all(InitialVectorType type_max = InitialVectorType::Alternating,
        InitialVectorType type_min = InitialVectorType::UnitConstant) {
        OutputHandler::print_header(ctx);
        compute_eigenvalues(type_max, type_min);
        auto solution = solve_seidel();
        OutputHandler::save_full_solution(ctx, "solution.bin", solution);
        OutputHandler::plot_full_solution(ctx);
    }

    /**
     * @brief Runs the Seidel method with saving and visualization of the solution.
     */
    [[maybe_unused]]
    void run_seidel() {
        OutputHandler::print_header(ctx);
        auto solution = solve_seidel();
        OutputHandler::save_full_solution(ctx, "solution.bin", solution);
        OutputHandler::plot_full_solution(ctx);
    }

    /**
     * @brief Runs the eigenvalue computation.
     * @param type_max Type of initial vector for the maximum eigenvalue.
     * @param type_min Type of initial vector for the minimum eigenvalue.
     */
    [[maybe_unused]]
    void run_eigen(InitialVectorType type_max = InitialVectorType::Alternating,
        InitialVectorType type_min = InitialVectorType::UnitConstant) {
        OutputHandler::print_header(ctx);
        compute_eigenvalues(type_max, type_min);
    }

private:
    /**
     * @brief Creates a Seidel solver and runs it.
     * @return Solution vector.
     */
    std::vector<double> solve_seidel() {
        std::cout << "\n--- Seidel method ---\n";
        SeidelSolver seidel(ctx.mesh, ctx.coeff_x, ctx.coeff_y);
        Clock clock;
        auto solution = seidel.solve(ctx.tol_seidel, ctx.max_iter_seidel);
        std::cout << "Time: " << clock.elapsed().count() << " s\n";
        return solution;
    }

    /**
     * @brief Generates an initial vector for eigenvalue iterations.
     * @param type Vector type.
     * @return Vector of size total_nodes() with zeros on the boundary and interior nodes filled.
     */
    std::vector<double> generate_initial_vector(InitialVectorType type) const {
        const std::size_t n = ctx.mesh.nodes_per_side();
        std::vector<double> v(ctx.mesh.total_nodes(), 0.0);

        switch (type) {
        case InitialVectorType::Random: {
            std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<double> dist(-1.0, 1.0);
            ctx.mesh.for_each_inner(
                [&](std::size_t, std::size_t, std::size_t idx) {
                    v[idx] = dist(rng);
                }
            );
            break;
        }
        case InitialVectorType::UnitConstant: {
            ctx.mesh.for_each_inner(
                [&](std::size_t, std::size_t, std::size_t idx) {
                    v[idx] = 1.0;
                }
            );
            break;
        }
        case InitialVectorType::Alternating: {
            for (std::size_t j = 1; j < n - 1; ++j) {
                for (std::size_t i = 1; i < n - 1; ++i) {
                    double val = ((i + j) % 2 == 0) ? 1.0 : -1.0;
                    v[ctx.mesh.full_index(i, j)] = val;
                }
            }
            break;
        }
        }
        return v;
    }

    /**
     * @brief Computes the maximum and minimum eigenvalues.
     * @param type_max Initial vector type for maximum.
     * @param type_min Initial vector type for minimum.
     * @details Prints results and theoretical values.
     */
    void compute_eigenvalues(InitialVectorType type_max, InitialVectorType type_min) {
        std::cout << "\n--- Eigenvalue computation ---\n";
        LaplacianOperator laplacian(ctx.mesh, ctx.coeff_x, ctx.coeff_y);

        std::vector<double> initial_max = generate_initial_vector(type_max);
        std::vector<double> initial_min = generate_initial_vector(type_min);

        Clock clock;

        PowerMethod::Params params(laplacian, ctx.mesh, initial_max,
            ctx.tol_power, ctx.max_iter_power);

        PowerMethod power_method;
        auto [lambda_max, iter_max] = power_method.compute_max(params);
        Clock::Elapsed elapsed_max = clock.restart();

        ShiftedPowerMethod shifted_method(lambda_max);
        params.x = initial_min;
        auto [lambda_min, iter_min] = shifted_method.compute_min(params);
        Clock::Elapsed elapsed_min = clock.elapsed();

        double h = ctx.mesh.step();
        auto [max_theor, min_theor] = eigen_theory::lambda_max_min(ctx.coeff_x, ctx.coeff_y, h);

        std::cout << "\n--- Max eigenvalue ---\n";
        std::cout << "Computed: " << lambda_max << "\n";
        std::cout << "Theoretical: " << max_theor << "\n";
        std::cout << "Difference: " << std::abs(lambda_max - max_theor) << "\n";
        std::cout << "Iterations: " << iter_max << "\n";
        std::cout << "Time: " << elapsed_max.count() << " s\n";

        std::cout << "\n--- Min eigenvalue ---\n";
        std::cout << "Computed: " << lambda_min << "\n";
        std::cout << "Theoretical: " << min_theor << "\n";
        std::cout << "Difference: " << std::abs(lambda_min - min_theor) << "\n";
        std::cout << "Iterations: " << iter_min << "\n";
        std::cout << "Time: " << elapsed_min.count() << " s\n";
    }

private:
    Context ctx;
};

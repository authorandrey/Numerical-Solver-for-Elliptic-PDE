#pragma once
#include "context.h"
#include <fstream>
#include <iostream>
#include <iomanip>

/**
 * @class OutputHandler
 * @brief Static class for outputting results: printing, saving to binary file, plotting.
 */
class OutputHandler {
public:
    /**
     * @brief Prints solver parameters to console.
     * @param ctx Context with parameters.
     */
	static void print_header(const Context& ctx) {
        std::cout << std::fixed << std::setprecision(12);
        std::cout << "========================================\n";
        std::cout << "Solving -a u_xx - b u_yy = f on [0,1]^2\n";
        std::cout << "a = " << ctx.a << ", b = " << ctx.b << "\n";
        std::cout << "h = " << ctx.h << ", h^2 = " << ctx.h * ctx.h << "\n";
        std::cout << "Tolerance (power) = " << ctx.tol_power << ", max iterations = " << ctx.max_iter_power << "\n";
        std::cout << "Tolerance (seidel) = " << ctx.tol_seidel << ", max iterations = " << ctx.max_iter_seidel << "\n";
        std::cout << "========================================\n";
    }

    /**
     * @brief Saves solution vector to a binary file.
     * @param ctx Context.
     * @param filename Output file name.
     * @param solution Solution vector (size total_nodes()).
     * @throws std::runtime_error if file cannot be opened.
     * @details Format: first two 8-byte integers (nx, ny), then nx*ny double values in row-major order.
     */
    static void save_full_solution(
        const Context& ctx,
        const std::string& filename,
        const std::vector<DataType> solution
    ) {
        std::size_t nx = ctx.mesh.nodes_per_side();
        std::size_t ny = ctx.mesh.nodes_per_side();
        std::ofstream out(filename, std::ios::binary);
        if (!out) {
            throw std::runtime_error("Error: cannot open file " + filename + " for writing.\n");
        }
        out.write(reinterpret_cast<const char*>(&nx), sizeof(nx));
        out.write(reinterpret_cast<const char*>(&ny), sizeof(ny));
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                double val = solution[ctx.mesh.full_index(i, j)];
                out.write(reinterpret_cast<const char*>(&val), sizeof(val));
            }
        }
        out.close();
        std::cout << "Full solution saved to " << filename << " ("
            << nx << "x" << ny << " points)\n";
    }

    /**
     * @brief Launches a Python script to visualize the solution.
     * @details Calls "python plot_full.py".
     * @throws std::runtime_error if system call fails.
     */
    static void plot_full_solution() {
        std::cout << "Launching Python script to plot full solution...\n";
        int ret = std::system("python plot_full.py");
        if (ret != 0) {
            throw std::runtime_error("Warning: Python script execution failed (code " + std::to_string(ret) + ").\n");
        }
    }
};
#pragma once
#include <stdexcept>
#include <functional>
#include "core/types.h"

/**
 * @class UniformMesh
 * @brief Represents a uniform rectangular grid on [0,1]x[0,1] with a given step.
 * @details The grid is square, nodes at points (i*h, j*h) where i,j = 0,...,N-1,
 *          h – actual step, N = 1/h + 1.
 */
class UniformMesh {
public:
    /**
     * @brief Constructor with desired step.
     * @param step Desired step. Must be in (0, 0.5]. The step
     *             is adjusted so that the number of intervals is integer: step = 1/(N-1).
     * @throws std::invalid_argument if step <= 0 or step > 0.5.
     */
    explicit UniformMesh(DataType step) : m_step(step) {
        if (step <= 0.0 || step > 0.5)
            throw std::invalid_argument("Step must be in (0, 0.5]");
        m_nodes_per_side = static_cast<std::size_t>(1.0 / step) + 1;
        m_step = 1.0 / (m_nodes_per_side - 1);
    }

    /**
     * @brief Returns the number of nodes per side.
     * @return Number of nodes N.
     */
    std::size_t nodes_per_side() const { return m_nodes_per_side; }

    /**
     * @brief Returns the grid step.
     * @return Step h.
     */
    DataType step() const { return m_step; }

    /**
     * @brief Computes the linear index of node (i,j).
     * @param i Index along x (0-based).
     * @param j Index along y (0-based).
     * @return Linear index = i + j * nodes_per_side().
     */
    std::size_t full_index(std::size_t i, std::size_t j) const {
        return i + j * m_nodes_per_side;
    }

    /**
     * @brief Checks whether a node is interior (not on the boundary).
     * @param i Index along x.
     * @param j Index along y.
     * @return true if 0 < i < N-1 and 0 < j < N-1.
     */
    bool is_inner(std::size_t i, std::size_t j) const {
        return i > 0 && i < m_nodes_per_side - 1 &&
            j > 0 && j < m_nodes_per_side - 1;
    }

    /**
     * @brief Returns the total number of interior nodes.
     * @return (N-2)*(N-2).
     */
    std::size_t inner_nodes() const {
        std::size_t n = m_nodes_per_side - 2;
        return n * n;
    }

    /**
     * @brief Returns the total number of nodes, including boundary.
     * @return N*N.
     */
    std::size_t total_nodes() const {
        return m_nodes_per_side * m_nodes_per_side;
    }

    /**
     * @brief Returns x-coordinate of node i.
     * @param i Index along x.
     * @return i * step.
     */
    DataType x(std::size_t i) const { return i * m_step; }
    /**
     * @brief Returns y-coordinate of node j.
     * @param j Index along y.
     * @return j * step.
     */
    DataType y(std::size_t j) const { return j * m_step; }

    using IndexIterationFunction = std::function<void(size_t, size_t, size_t)>;

    /**
     * @brief Applies a function to all interior nodes.
     * @param f Function taking (i, j, linear_index).
     * @note The function is called for each interior node in row-major order.
     */
    void for_each_inner(const IndexIterationFunction& f) const {
        for (std::size_t j = 1; j < m_nodes_per_side - 1; ++j)
            for (std::size_t i = 1; i < m_nodes_per_side - 1; ++i)
                f(i, j, full_index(i, j));
    }

private:
    DataType m_step;
    std::size_t m_nodes_per_side;
};

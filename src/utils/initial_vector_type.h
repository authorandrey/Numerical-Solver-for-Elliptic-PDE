#pragma once

/**
 * @enum InitialVectorType
 * @brief Type of initial vector for eigenvalue iterations.
 */
enum class InitialVectorType {
    Random,         ///< Random values in [-1,1].
    UnitConstant,   ///< All interior nodes equal to 1.
    Alternating     ///< Alternating sign (-1)^(i+j).
};

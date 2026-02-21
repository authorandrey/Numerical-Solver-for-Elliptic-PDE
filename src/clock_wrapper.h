#pragma once
#include "context.h"
#include <chrono>

/**
 * @class Clock
 * @brief Simple RAII timer for measuring elapsed time.
 * @details Uses std::chrono::high_resolution_clock.
 */
class Clock {
public:
    using TimePoint = std::chrono::steady_clock::time_point;
    using ChronoClock = std::chrono::high_resolution_clock;
    using Elapsed = std::chrono::duration<DataType>;

    /**
     * @brief Constructor starts the timer.
     */
    Clock() : m_start(ChronoClock::now()) { }

    /**
     * @brief Returns the time elapsed since construction or the last restart.
     * @return Duration in seconds (DataType).
     */
    Elapsed elapsed() {
        auto end = ChronoClock::now();
        return end - m_start;
    }

    /**
     * @brief Restarts the timer and returns the elapsed time.
     * @return Duration in seconds (DataType).
     */
    Elapsed restart() {
        Elapsed elapsed = std::move(this->elapsed());
        m_start = ChronoClock::now();
        return elapsed;
    }

private:
    TimePoint m_start;
};
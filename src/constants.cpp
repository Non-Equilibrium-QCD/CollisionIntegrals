#pragma once


#include <cstddef>
constexpr double g = 1.0;
constexpr double dA = 8.0;
constexpr double dF = 3.0;
constexpr double CA = 3.0;
constexpr double CF = 4.0 / 3.0;
constexpr double mDSqr = 1.0;
constexpr double mQSqr = 1.0;
constexpr double nuG = 2.0 * dA;
constexpr double nuQ = 2.0 * dF;

template <typename T, size_t exp>
constexpr T power_recursive(T base) {

    static_assert(exp >= 0, "Exponent must be non-negative");

    // Recursive step: x^n = x * x^(n-1)
    T result = 1;
    for (int i = 0; i < exp; ++i) {
        result *= base;
    }
    return result;
}

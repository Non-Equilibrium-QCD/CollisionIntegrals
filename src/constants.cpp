#pragma once


#include <cmath>
#include <cstddef>
#include <string>
constexpr double dA = 8.0;
constexpr double dF = 3.0;
constexpr double CA = 3.0;
constexpr double CF = 4.0 / 3.0;
constexpr double mDSqr = 1.0;
constexpr double mQSqr = 1.0;
constexpr double nuG = 2.0 * dA;
constexpr double nuQ = 2.0 * dF;
constexpr double Nc = 3.0;
constexpr double Nf = 3.0;
constexpr double lambda = 10.0;
const double g = sqrt(lambda / Nc);

template <typename T, size_t exp>
constexpr T power_recursive(T base) {

    static_assert(exp >= 0, "Exponent must be non-negative");

    // Recursive step: x^n = x * x^(n-1)
    T result = 1;
    for (size_t i = 0; i < exp; ++i) {
        result *= base;
    }
    return result;
}


std::string ReadFileName(int argc, char *argv[], const std::string &defaultName) {
    if (argc >= 2) {
        return std::string(argv[1]);
    } else {
        return defaultName;
    }
}

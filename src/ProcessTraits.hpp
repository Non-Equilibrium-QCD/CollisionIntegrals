#pragma once

/**
 * @file ProcessTraits.hpp
 * @brief Compile-time process definitions for collision integrals.
 *
 * Each process defines:
 * - Distribution functions for particles 1-4
 * - Statistical factor
 * - Matrix element squared
 */

// ============================================================================
// Process Traits
// ============================================================================

/**
 * @brief Process trait template.
 *
 * Specialize this for each scattering process to define:
 * - dist1, dist2, dist3, dist4: distribution functions
 * - stat: statistical factor
 * - matrix: matrix element squared
 */
template<typename ProcessTag>
struct ProcessTraits;

// ============================================================================
// Process Tags Reference Implementation
// ============================================================================

/// Tag for xx -> xx process
struct x1x2_x3x4 {};

namespace Distribution {
// static should make the function be internal to this file
inline static double nop(double p, double cosTheta, double phi) {
    (void)p; (void)cosTheta; (void)phi;
    return 0.0;
}
}

namespace MatrixElement {
inline static double nop(double s, double t, double u,
                    double qt, double qu,
                    double mDSqr, double mQSqr) {
    (void)s; (void)t; (void)u;
    (void)qt; (void)qu;
    (void)mDSqr; (void)mQSqr;
    return 0.0;
}
} // namespace MatrixElement

namespace StatisticalFactor {
inline static double nop(double f1, double f2, double f3, double f4) {
    (void)f1; (void)f2; (void)f3; (void)f4;
    return 0.0;
}
} // namespace StatisticalFactor

template<>
struct ProcessTraits<x1x2_x3x4> {
    static constexpr auto dist1 = Distribution::nop;
    static constexpr auto dist2 = Distribution::nop;
    static constexpr auto dist3 = Distribution::nop;
    static constexpr auto dist4 = Distribution::nop;
    static constexpr auto stat = StatisticalFactor::nop;
    static constexpr auto matrix = MatrixElement::nop;
    static constexpr double nuA = 1.0;
};

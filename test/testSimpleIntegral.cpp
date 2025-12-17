#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../third_party/doctest/doctest.h"

#include "../src/constants.cpp"
#include "../src/Process.hpp"
#include "../src/GSLVEGAS.cpp"
#include <cmath>
#include <fmt/core.h>

// ============================================================================
// SIMPLE TEST INTEGRAL
// ============================================================================
//
// From Scalar.md: A simplified collision integral with:
//   - Constant matrix element: |M|² = 1
//   - Statistical factor: F[f] = f₂ = exp(-p₂²)
//   - No angular dependence (integrates to (2π)² for Φ₁q and Φ₂q)
//
// The integral is:
//   C[f₁] = 1/(2ν₁) × 1/(2π)⁴ × 1/(16p₁²) × ∫dp₂ ∫dq ∫dω/(2π) × exp(-p₂²)
//
// Analytical result (from Sage):
//   C[f₁] = p₁ / (2ν₁ × (2π)⁵ × 16p₁²) = 1 / (2ν₁ × (2π)⁵ × 16p₁)
//
// ============================================================================

namespace SimpleIntegral {

// Simple statistical factor: just f2
inline double statSimple(double f1, double f2, double f3, double f4) {
    (void)f1; (void)f3; (void)f4;  // unused
    return f2;
}

// Constant matrix element
inline double matrixConstant(double s, double t, double u,
                             double qt, double qu,
                             double mDSqr, double mQSqr) {
    (void)s; (void)t; (void)u; (void)qt; (void)qu;
    (void)mDSqr; (void)mQSqr;
    return 1.0;
}

// Gaussian distribution for p2
inline double gaussianDist(double p, double cosTheta, double phi) {
    (void)cosTheta; (void)phi;
    return std::exp(-p * p);
}

// Unit distribution (returns 1)
inline double unitDist(double p, double cosTheta, double phi) {
    (void)p; (void)cosTheta; (void)phi;
    return 1.0;
}

}  // namespace SimpleIntegral

// Process tag for the simple test integral
struct simple_test {};

template<>
struct ProcessTraits<simple_test> {
    static constexpr auto dist1 = SimpleIntegral::unitDist;
    static constexpr auto dist2 = SimpleIntegral::gaussianDist;
    static constexpr auto dist3 = SimpleIntegral::unitDist;
    static constexpr auto dist4 = SimpleIntegral::unitDist;
    static constexpr auto stat = SimpleIntegral::statSimple;
    static constexpr auto matrix = SimpleIntegral::matrixConstant;
    static constexpr double nuA = 1.0;  // Simplified degeneracy
};

// Include Integral.cpp after defining the process traits
#include "../src/Integral.cpp"

// ============================================================================
// Analytical Result
// ============================================================================
//
// C[f_1] = 1 / (2 × nuA × (2π)^3 × 16 × p₁)
//
inline double analyticalResult(double p1, double nuA = 1.0) {
    const double twoPi = power_recursive<double, 3>(2.0 * M_PI);
    return 1.0 / (2.0 * nuA * twoPi * 16.0 * p1);
}

// ============================================================================
// TEST SUITE
// ============================================================================

TEST_SUITE("Simple Test Integral") {

    TEST_CASE("VEGAS result matches analytical value within tolerance") {
        double p1_values[] = {0.5, 1.0, 2.0};

        for (double p1 : p1_values) {
            GSLARGS args;
            args.p1 = p1;
            args.cosTheta1 = 0.0;
            args.phi1 = 0.0;
            args.fct = CollisionIntegralQCD::CollisionIntegral<simple_test>;

            GSLVEGAS vegas(5);
            double result, error;

            // Use more calls for better precision
            vegas.integrate(args, 50000, &result, &error);

            double expected = analyticalResult(p1);

            MESSAGE("p1 = ", p1);
            MESSAGE("  VEGAS:      ", result, " +/- ", error);
            MESSAGE("  Analytical: ", expected);
            MESSAGE("  Difference: ", std::abs(result - expected));
            MESSAGE("  Ratio:      ", result / expected);

            // Check within 3 sigma or 10% relative error (whichever is larger)
            double tolerance = std::max(3.0 * error, 0.1 * std::abs(expected));
            CHECK(std::abs(result - expected) < tolerance);
        }
    }

}

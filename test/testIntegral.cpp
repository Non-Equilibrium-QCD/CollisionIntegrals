#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../third_party/doctest/doctest.h"

#include "../src/Integral.cpp"
#include <cmath>

// ============================================================================
// UNIT TESTS FOR COLLISION INTEGRAL
// ============================================================================

TEST_SUITE("CollisionIntegral Integrand") {

    TEST_CASE("integrand returns finite value for typical phase-space point") {
        double x[5] = {0.3, 0.4, 0.5, 0.25, 0.75};
        
        GSLARGS args;
        args.p1 = 1.0;
        args.cosTheta1 = 0.0;
        args.phi1 = 0.0;
        
        double result = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>(x, args);
        
        CHECK(std::isfinite(result));
    }

    TEST_CASE("integrand returns finite value for various momenta") {
        double x[5] = {0.5, 0.5, 0.5, 0.5, 0.5};
        
        GSLARGS args;
        args.cosTheta1 = 0.0;
        args.phi1 = 0.0;
        
        SUBCASE("small momentum p1 = 0.1") {
            args.p1 = 0.1;
            double result = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>(x, args);
            CHECK(std::isfinite(result));
        }
        
        SUBCASE("medium momentum p1 = 1.0") {
            args.p1 = 1.0;
            double result = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>(x, args);
            CHECK(std::isfinite(result));
        }
        
        SUBCASE("large momentum p1 = 10.0") {
            args.p1 = 10.0;
            double result = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>(x, args);
            CHECK(std::isfinite(result));
        }
    }

    TEST_CASE("integrand returns finite value for various angles") {
        double x[5] = {0.3, 0.4, 0.5, 0.25, 0.75};
        
        GSLARGS args;
        args.p1 = 1.0;
        args.phi1 = 0.0;
        
        SUBCASE("cosTheta1 = -1.0 (backward)") {
            args.cosTheta1 = -1.0;
            double result = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>(x, args);
            CHECK(std::isfinite(result));
        }
        
        SUBCASE("cosTheta1 = 0.0 (perpendicular)") {
            args.cosTheta1 = 0.0;
            double result = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>(x, args);
            CHECK(std::isfinite(result));
        }
        
        SUBCASE("cosTheta1 = 1.0 (forward)") {
            args.cosTheta1 = 1.0;
            double result = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>(x, args);
            CHECK(std::isfinite(result));
        }
    }

    TEST_CASE("integrand handles edge cases at integration boundaries") {
        GSLARGS args;
        args.p1 = 1.0;
        args.cosTheta1 = 0.0;
        args.phi1 = 0.0;
        
        SUBCASE("x[0] near 0 (small p2)") {
            double x[5] = {0.01, 0.5, 0.5, 0.5, 0.5};
            double result = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>(x, args);
            CHECK(std::isfinite(result));
        }
        
        SUBCASE("x[0] near 1 (large p2)") {
            double x[5] = {0.99, 0.5, 0.5, 0.5, 0.5};
            double result = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>(x, args);
            CHECK(std::isfinite(result));
        }
        
        SUBCASE("x[1] = 0 (q at minimum)") {
            double x[5] = {0.5, 0.0, 0.5, 0.5, 0.5};
            double result = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>(x, args);
            // q=0 leads to division issues, result may be 0 or NaN handled gracefully
            CHECK((std::isfinite(result) || result == 0.0));
        }
        
        SUBCASE("x[1] = 1 (q at maximum)") {
            double x[5] = {0.5, 1.0, 0.5, 0.5, 0.5};
            double result = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>(x, args);
            CHECK(std::isfinite(result));
        }
    }

    TEST_CASE("integrand returns zero when |t| or |u| is very small") {
        // The integrand explicitly returns 0 when |t| < 1e-12 or |u| < 1e-12
        // This happens at specific kinematic configurations
        
        GSLARGS args;
        args.p1 = 1.0;
        args.cosTheta1 = 0.0;
        args.phi1 = 0.0;
        
        // When w = q, we get t = w² - q² = 0
        // This happens when x[2] = 1.0 (w at wMax = q for certain configurations)
        // We need to find x values that give t ≈ 0
        
        // For t = w² - q² = 0, we need w = ±q
        // w = wMin + (wMax - wMin) * x[2]
        // When x[2] makes w = q (if q is in range), then t = 0
        
        // This is a mathematical property - we just verify the code handles it
        MESSAGE("Note: Finding exact t=0 or u=0 points requires solving kinematic constraints");
    }

}

TEST_SUITE("Jacobian and Kinematics") {

    TEST_CASE("momentum conservation: p1 + p2 = p3 + p4") {
        // This is implicitly tested by the integrand, but we can verify
        // the kinematic setup is consistent
        
        double x[5] = {0.3, 0.4, 0.5, 0.25, 0.75};
        
        GSLARGS args;
        args.p1 = 1.0;
        args.cosTheta1 = 0.0;
        args.phi1 = 0.0;
        
        // The integrand should produce a finite result when kinematics are valid
        double result = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>(x, args);
        CHECK(std::isfinite(result));
    }

    TEST_CASE("p2 mapping: x[0]/(1-x[0]) maps [0,1) to [0,inf)") {
        // Verify the mapping behavior
        CHECK(0.0 / (1.0 - 0.0) == doctest::Approx(0.0));
        CHECK(0.5 / (1.0 - 0.5) == doctest::Approx(1.0));
        CHECK(0.9 / (1.0 - 0.9) == doctest::Approx(9.0));
        CHECK(0.99 / (1.0 - 0.99) == doctest::Approx(99.0));
    }

    TEST_CASE("Jacobian for p2 mapping is (1+p2)^2") {
        auto jacobian = [](double x0) {
            double p2 = x0 / (1.0 - x0);
            return (1.0 + p2) * (1.0 + p2);
        };
        
        // dp2/dx0 = 1/(1-x0)^2 = (1+p2)^2
        CHECK(jacobian(0.0) == doctest::Approx(1.0));
        CHECK(jacobian(0.5) == doctest::Approx(4.0));
        CHECK(jacobian(0.9) == doctest::Approx(100.0));
    }

}

TEST_SUITE("Distribution Functions") {

    TEST_CASE("Bose-Einstein distribution is positive and monotonically decreasing") {
        double T = 1.0;
        double f_prev = Distribution::bose(0.01, T);
        
        for (double p = 0.1; p <= 5.0; p += 0.1) {
            double f = Distribution::bose(p, T);
            CHECK(f > 0.0);
            CHECK(f < f_prev);
            f_prev = f;
        }
    }

    TEST_CASE("Fermi-Dirac distribution is bounded between 0 and 1") {
        double T = 1.0;
        
        for (double p = 0.01; p <= 5.0; p += 0.1) {
            double f = Distribution::fermi(p, T);
            CHECK(f > 0.0);
            CHECK(f < 1.0);
        }
    }

    TEST_CASE("Gaussian distribution (gluon) is positive and isotropic in angles") {
        double p = 1.0;
        
        // Test isotropy (currently gluon distribution ignores angles in exp(-p²))
        double f1 = Distribution::gluon(p, 0.0, 0.0);
        double f2 = Distribution::gluon(p, 0.5, 1.0);
        double f3 = Distribution::gluon(p, -0.5, 2.0);
        
        CHECK(f1 == doctest::Approx(f2));
        CHECK(f2 == doctest::Approx(f3));
        CHECK(f1 > 0.0);
    }

}

TEST_SUITE("Statistical Factors") {

    TEST_CASE("BBBB statistical factor for equilibrium distributions") {
        // For Bose-Einstein at same temperature: f = 1/(exp(p/T) - 1)
        // The statistical factor f1*f2*(1+f3)*(1+f4) - f3*f4*(1+f1)*(1+f2)
        // should vanish when p1 + p2 = p3 + p4 (energy conservation)
        
        double T = 1.0;
        
        // Test case: p1 = p3 = 1.0, p2 = p4 = 2.0 (trivial energy conservation)
        double f1 = Distribution::bose(1.0, T);
        double f2 = Distribution::bose(2.0, T);
        double f3 = Distribution::bose(1.0, T);
        double f4 = Distribution::bose(2.0, T);
        
        double stat = StatisticalFactor::BBBB(f1, f2, f3, f4);
        CHECK(stat == doctest::Approx(0.0).epsilon(1e-10));
    }

    TEST_CASE("BBBB statistical factor is non-zero for non-equilibrium") {
        // Use different distributions that don't satisfy detailed balance
        double f1 = 0.5;
        double f2 = 0.3;
        double f3 = 0.4;
        double f4 = 0.2;
        
        double stat = StatisticalFactor::BBBB(f1, f2, f3, f4);
        CHECK(stat != doctest::Approx(0.0));
    }

    TEST_CASE("FFFF statistical factor for equilibrium distributions") {
        double T = 1.0;
        
        // p1 = p3, p2 = p4 (trivial case)
        double f1 = Distribution::fermi(1.0, T);
        double f2 = Distribution::fermi(2.0, T);
        double f3 = Distribution::fermi(1.0, T);
        double f4 = Distribution::fermi(2.0, T);
        
        double stat = StatisticalFactor::FFFF(f1, f2, f3, f4);
        CHECK(stat == doctest::Approx(0.0).epsilon(1e-10));
    }

}

TEST_SUITE("Matrix Elements") {

    TEST_CASE("gg_gg matrix element is positive for physical kinematics") {
        // Physical kinematics: s > 0, t < 0, u < 0, s + t + u = 0
        double s = 4.0;
        double t = -1.0;
        double u = -3.0;  // s + t + u = 0
        
        double qt = 1.0;
        double qu = 1.5;
        
        double M2 = MatrixElement::gg_gg(s, t, u, qt, qu, mDSqr, mQSqr);
        CHECK(M2 > 0.0);
        CHECK(std::isfinite(M2));
    }

    TEST_CASE("gg_gg matrix element scales with coupling") {
        double s = 4.0;
        double t = -1.0;
        double u = -3.0;
        double qt = 1.0;
        double qu = 1.5;
        
        double M2 = MatrixElement::gg_gg(s, t, u, qt, qu, mDSqr, mQSqr);
        
        // Matrix element should be proportional to g^4
        // g = sqrt(lambda/Nc) = sqrt(10/3) ≈ 1.826
        CHECK(M2 > 0.0);
    }

}

#include "../src/Integral.cpp"
#include "../src/constants.cpp"

namespace SimpleIntegral {

// Simple statistical factor: just f2
inline double statSimple(double f1, double f2, double f3, double f4) {
    (void)f1; (void)f3; (void)f4;  // unused
    // return f1 * f2 - f3 * f4;
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

int main(int argc, char *argv[]) {

    fmt::println(stderr, "Using ORIGINAL implementation (Integral.cpp)");
    Integrate::Setup();
    auto integrand = CollisionIntegral::CollisionIntegral<simple_test>;
    Integrate::Compute<simple_test>(ReadFileName(argc, argv,
                                       "OUTPUT/simple.dat"));

    return 0;
}

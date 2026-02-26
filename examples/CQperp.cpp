/**
 * @file CQperp.cpp
 * @brief Example: Compute C(q_perp) for arbitrary distribution and matrix element.
 *
 * This example computes the transverse momentum broadening kernel C(q_perp)
 * following the formalism in docs/source/CqPerp.md.
 *
 * Usage:
 *   make example NAME=CQperp
 *   ./examples/CQperp [output_file]
 */

#include "../src/CQper.cpp"

// ============================================================================
// Distribution Functions
// ============================================================================

namespace DistributionCQperp {

/// Bose-Einstein distribution for gluons
inline double bose(double p, double cosTheta, double phi) {
    (void)cosTheta; (void)phi;
    return 1.0 / std::expm1(p);
}

/// Fermi-Dirac distribution for quarks
inline double fermi(double p, double cosTheta, double phi) {
    (void)cosTheta; (void)phi;
    return 1.0 / (std::exp(p) + 1.0);
}

/// Non-equilibrium Gaussian distribution
inline double gaussian(double p, double cosTheta, double phi) {
    (void)cosTheta; (void)phi;
    return std::exp(-p * p);
}

// DISTRIBUTION fInit(p,Theta) //
double fGInit(double p, double CosTheta, double phi) {
    double Xi0 = 10.0;
    (void)phi;

    double A = 5.34; double Q0 = 1.8; double lambda = g * g * Nc;

    return 2.0 * A / lambda * (Q0 / p) *
           std::exp(-2.0 / 3.0 * (p * p) / (Q0 * Q0) *
                    (1.0 + (Xi0 * Xi0 - 1.0) * CosTheta * CosTheta))
           / std::sqrt((1.0 + (Xi0 * Xi0 - 1.0) * CosTheta * CosTheta));
}


} // namespace DistributionCQperp

// ============================================================================
// Matrix Elements for C(q_perp)
//
// Signature: double matrix(double s, double t,
//                          double qt, double mDSqr, double mQSqr)
// where s = 2*(p2 - p2_parallel), t = -q_perp^2, u = -s - t
// ============================================================================

namespace MatrixElementCQperp {

// SCREENING PARAMETERS //
static const double xi0g = std::exp(5.0 / 6.0) / (2.0 * M_SQRT2);
static const double xi0q = std::exp(1.0) / M_SQRT2;

/**
 * @brief Leading-Log matrix element for C(q_perp).
 *
 * In the p1 -> infinity limit with s = 2*(p2 - p2_parallel):
 * |M|^2_LL includes screening via tBar
 */
inline double LeadingLog(double s, double t, double w,
                         double qt,
                         double mDSqr, double mQSqr, double c1, double c2) {
    (void)mQSqr;
    (void)w;
    (void)c1;
    (void)c2;
    double tBar = t * (qt * qt + xi0g * xi0g * mDSqr) / (qt * qt);
    // s = 2*(p2 - p2_parallel), so s^2/4 = (p2 - p2_parallel)^2
    return 4.0 * g * g * g * g * dA * CA * CA * (
               s * s / (4.0 * tBar * tBar)
           );
}

/**
 * @brief Leading-Log matrix element for C(q_perp).
 *
 * In the p1 -> infinity limit with s = 2*(p2 - p2_parallel):
 * |M|^2_LL includes screening via tBar
 */
inline double HTL(double s, double t, double w,
                  double q,
                  double mDSqr, double mQSqr, double c1, double c2) {
    (void)s;
    (void)t;
    (void)mQSqr;
    double AHTL = q * q
                  + mDSqr * (1 + (w / (2 * q)) *
                             std::log((q - w) / (q + w)));
    double BHTL = - (mDSqr * w / (2 * q)) * M_PI;
    double CHTL = q * q - w * w + (mDSqr / 2) *
                  ((w * w / q * q) + ((w * w / q * q) - 1)
                   * (w / (2 * q)) * std::log((q - w) / (q + w)));
    double DHTL = (M_PI * mDSqr * w / (4 * q)) * (1 - (w * w / q * q));

    return c1 * c1 / (AHTL * AHTL + BHTL * BHTL) + c2 * c2 /
           (CHTL * CHTL + DHTL * DHTL) -
           (2 * c1 * c2 * (AHTL * CHTL + BHTL * DHTL)) / ((AHTL * AHTL + BHTL * BHTL) *
                   (CHTL * CHTL + DHTL * DHTL));

}

} // namespace MatrixElementCQperp

// ============================================================================
// Statistical Factors
// ============================================================================

namespace StatisticalCQperp {

/// Bose enhancement: C_A * f2 * (1 + f4)
inline double BoseEnhancement(double f1, double f2, double f3, double f4) {
    (void)f1; (void)f3;
    return CA * f2 * (1.0 + f4);
}

/// Fermi blocking: (N_f/2) * f2 * (1 - f4)
inline double FermiBlocking(double f1, double f2, double f3, double f4) {
    (void)f1; (void)f3;
    return (Nf / 2.0) * f2 * (1.0 - f4);
}

/// Full QCD statistical factor (gluons only for simplicity)
inline double FullQCD(double f1, double f2, double f3, double f4) {
    (void)f1; (void)f3;
    // C_A * f_g(p2) * (1 + f_g(p4))
    // For full QCD, would add quark terms with separate distributions
    return CA * f2 * (1.0 + f4);
}

} // namespace StatisticalCQperp

// ============================================================================
// Process Definitions
// ============================================================================

/// Tag for gluon-gluon C(q_perp) with thermal distribution (Leading-Log)
struct CQperp_thermal_LL {};

template<>
struct ProcessTraits<CQperp_thermal_LL> {
    static constexpr auto dist1 =
        DistributionCQperp::bose;  // not used in C(q_perp)
    static constexpr auto dist2 = DistributionCQperp::bose;  // medium particle
    static constexpr auto dist3 = DistributionCQperp::bose;  // not used
    static constexpr auto dist4 =
        DistributionCQperp::bose;  // scattered medium particle
    static constexpr auto stat  = StatisticalCQperp::BoseEnhancement;
    static constexpr auto matrix = MatrixElementCQperp::LeadingLog;
    static constexpr double nuA = nuG;
};

/// Tag for non-equilibrium Gaussian distribution
struct CQperp_gaussian_LL {};

template<>
struct ProcessTraits<CQperp_gaussian_LL> {
    static constexpr auto dist1 = DistributionCQperp::gaussian;
    static constexpr auto dist2 = DistributionCQperp::gaussian;
    static constexpr auto dist3 = DistributionCQperp::gaussian;
    static constexpr auto dist4 = DistributionCQperp::gaussian;
    static constexpr auto stat  = StatisticalCQperp::BoseEnhancement;
    static constexpr auto matrix = MatrixElementCQperp::LeadingLog;
    static constexpr double nuA = nuG;
};

/// Tag for anisotropic distribution
struct CQperp_anisotropic_LL {};

template<>
struct ProcessTraits<CQperp_anisotropic_LL> {
    static constexpr auto dist1 = DistributionCQperp::fGInit;
    static constexpr auto dist2 = DistributionCQperp::fGInit;
    static constexpr auto dist3 = DistributionCQperp::fGInit;
    static constexpr auto dist4 = DistributionCQperp::fGInit;
    static constexpr auto stat  = StatisticalCQperp::BoseEnhancement;
    // static constexpr auto matrix = MatrixElementCQperp::LeadingLog;
    static constexpr auto matrix = MatrixElementCQperp::HTL;
    static constexpr double nuA = nuG;
};

// ============================================================================
// Main
// ============================================================================

int main(int argc, char *argv[]) {
    // Parse command line arguments
    std::string outputFile = ReadFileName(argc, argv,
                                          "OUTPUT/CQperp_thermal_LL.dat");

    fmt::println(stderr, "Computing C(q_perp) with Leading-Log matrix element");
    fmt::println(stderr, "Output file: {}", outputFile);
    fmt::println(stderr, "Using {} OpenMP threads", omp_get_max_threads());

    // Initialize integrators
    IntegrateCQperp::Setup();

    // Compute C(q_perp) for thermal distribution with Leading-Log
    // Arguments: OutputFile, p1 (momentum), cosTheta1 (angle)
    // IntegrateCQperp::Compute<CQperp_thermal_LL>(outputFile, 1.0, 1.0);
    IntegrateCQperp::Compute<CQperp_anisotropic_LL>(outputFile, 1.0, 0.0);

    fmt::println(stderr, "Done!");
    return EXIT_SUCCESS;
}

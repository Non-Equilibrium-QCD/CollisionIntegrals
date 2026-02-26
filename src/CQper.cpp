#pragma once

#include <cmath>
#include "constants.cpp"
#include "ProcessTraits.hpp"
#include <gsl/gsl_monte_vegas.h>
#include <vector>
#include "GSLVEGAS.cpp"
#include <omp.h>
#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/color.h>
// #define CQPERP_LINEAR_P2 // Uncomment to use linear sampling for p2 instead of transform sampling

namespace CQperpIntegral {
constexpr size_t dimensions = 3; // p2, omega, phi2q

struct GSLARGS {
    double qPerp;           // transverse momentum transfer
    double Phi1Q;           // azimuthal angle of q relative to p1
    double p1, cosTheta1, phi1; // particle 1 momentum
    double (*fct)(double *, const GSLARGS &);
};

/**
 * @brief Computes C(q_perp) for 2->2 processes using Monte Carlo integration.
 *
 * This implements the integral from CqPerp.md:
 * C(q_perp) = (1/(2pi)^2) * (1/q_perp^4) * integral[dw dp2 dphi2q] |M|^2 * F[f]
 *
 * The integration variables are:
 * - x[0]: maps to p2 in [(q-w)/2, infinity)
 * - x[1]: maps to omega in (-infinity, +infinity) via tan transform
 * - x[2]: maps to phi2q in [0, 2*pi)
 *
 * @tparam ProcessTag Tag type identifying the process.
 *         ProcessTraits<ProcessTag> must define:
 *         - dist2, dist4: distribution functions for particles 2 and 4
 *         - stat: statistical factor function F[f](f1, f2, f3, f4)
 *         - matrix: matrix element function |M|^2(s, t, u, qt, qu, mDSqr, mQSqr)
 * @param x Array of random variables for sampling.
 * @param args Struct containing qPerp and p1 angular coordinates.
 * @return The value of the integrand for C(q_perp).
 */
template<typename ProcessTag>
inline double CQperpIntegrand(double *x, const GSLARGS& args) {
    using Traits = ProcessTraits<ProcessTag>;

    (void)args.p1; // p1 -> infinity limit
    double cosTheta1 = args.cosTheta1;
    double Sin1 = std::sqrt(1.0 - cosTheta1 * cosTheta1);
    double phi1 = args.phi1;
    // Azimuthal angle of q relative to p1 (from args)
    double qPerp = args.qPerp;
    double Phi1Q = args.Phi1Q;

    // SET BASIS VECTORS WRT p1 //
    double ep1[3], ep2[3], ep3[3];

    ep1[0] = Sin1 * std::cos(phi1);
    ep1[1] = Sin1 * std::sin(phi1);
    ep1[2] = cosTheta1;

    ep2[0] = cosTheta1 * std::cos(phi1);
    ep2[1] = cosTheta1 * std::sin(phi1);
    ep2[2] = -Sin1;

    ep3[0] = -std::sin(phi1);
    ep3[1] = std::cos(phi1);
    ep3[2] = 0.0;

    // SAMPLE omega from -infinity to +infinity via tan transform //
    // x[1] in (0,1) -> omega in (-inf, +inf)
    double omega = std::tan(M_PI * (x[1] - 0.5));
    double dOmega_dx = M_PI / (std::cos(M_PI * (x[1] - 0.5)) * std::cos(M_PI *
                               (x[1] - 0.5)));

    // Compute q = sqrt(qPerp^2 + omega^2)
    double q = std::sqrt(qPerp * qPerp + omega * omega);

    if (q < 1e-12) {
        return 0.0;
    }

    // SAMPLE p2 from (q - omega)/2 to infinity (or p2Max) //
    double p2Min = std::max(0.0, (q - omega) / 2.0);
    double p2, dp2_dx;

    #ifdef CQPERP_LINEAR_P2
    // Linear sampling: p2 in [p2Min, p2Max]
    constexpr double p2Max = 100.0;
    p2 = p2Min + (p2Max - p2Min) * x[0];
    dp2_dx = (p2Max - p2Min);
    #else

    // Transform sampling: p2 = p2Min + y/(1-y) for y in [0,1) -> p2 in [p2Min, inf)
    if (x[0] >= 1.0 - 1e-12) {
        return 0.0;
    }

    double y = x[0];
    p2 = p2Min + y / (1.0 - y);
    dp2_dx = 1.0 / ((1.0 - y) * (1.0 - y));
    #endif

    // SAMPLE Phi2q //
    double Phi2Q = 2.0 * M_PI * x[2];
    double dPhi2Q = 2.0 * M_PI;

    // Compute Cos2Q from delta function constraint:
    // cos(theta_2q) = omega/q + (omega^2 - q^2)/(2*p2*q)
    //               = omega/q - qPerp^2/(2*p2*q)
    double Cos2Q = omega / q - qPerp * qPerp / (2.0 * p2 * q);

    // Check if Cos2Q is in valid range [-1, 1]
    if (std::abs(Cos2Q) > 1.0) {
        return 0.0;
    }

    double Sin2Q = std::sqrt(1.0 - Cos2Q * Cos2Q);

    // Compute Cos1Q = omega/q (angle between p1 and q)
    double Cos1Q = omega / q;
    double Sin1Q = std::sqrt(1.0 - Cos1Q * Cos1Q);


    // SET q VECTOR using basis vectors wrt p1 //
    double qx = q * (Cos1Q * ep1[0] + Sin1Q * std::cos(Phi1Q) * ep2[0] + Sin1Q *
                     std::sin(Phi1Q) * ep3[0]);
    double qy = q * (Cos1Q * ep1[1] + Sin1Q * std::cos(Phi1Q) * ep2[1] + Sin1Q *
                     std::sin(Phi1Q) * ep3[1]);
    double qz = q * (Cos1Q * ep1[2] + Sin1Q * std::cos(Phi1Q) * ep2[2] + Sin1Q *
                     std::sin(Phi1Q) * ep3[2]);

    // DETERMINE ANGLES of q vector //
    double CosQ = qz / q;
    double SinQ = std::sqrt(1.0 - CosQ * CosQ);
    double PhiQ = std::atan2(qy, qx);

    // SET BASIS VECTORS WRT q //
    double eq1[3], eq2[3], eq3[3];

    eq1[0] = SinQ * std::cos(PhiQ);
    eq1[1] = SinQ * std::sin(PhiQ);
    eq1[2] = CosQ;

    eq2[0] = CosQ * std::cos(PhiQ);
    eq2[1] = CosQ * std::sin(PhiQ);
    eq2[2] = -SinQ;

    eq3[0] = -std::sin(PhiQ);
    eq3[1] = std::cos(PhiQ);
    eq3[2] = 0.0;

    // SET p2 VECTOR //
    double p2x = p2 * (Cos2Q * eq1[0] + Sin2Q * std::cos(Phi2Q) * eq2[0] + Sin2Q *
                       std::sin(Phi2Q) * eq3[0]);
    double p2y = p2 * (Cos2Q * eq1[1] + Sin2Q * std::cos(Phi2Q) * eq2[1] + Sin2Q *
                       std::sin(Phi2Q) * eq3[1]);
    double p2z = p2 * (Cos2Q * eq1[2] + Sin2Q * std::cos(Phi2Q) * eq2[2] + Sin2Q *
                       std::sin(Phi2Q) * eq3[2]);

    // SET p4 VECTOR = p2 + q //
    double p4x = p2x + qx;
    double p4y = p2y + qy;
    double p4z = p2z + qz;
    double p4 = std::sqrt(p4x * p4x + p4y * p4y + p4z * p4z);

    // Compute p2_parallel = p2 . e^1_p1 (projection along p1 direction)
    double p2Parallel = p2x * ep1[0] + p2y * ep1[1] + p2z * ep1[2];

    // Compute Mandelstam variables in p1 -> infinity limit
    double s = 2.0 * (p2 - p2Parallel);
    double t = omega * omega - q * q; // = -qPerp^2

    // Momentum transfers
    double qt = q;  // |p1 - p3| = |q|

    // Matrix element from ProcessTraits
    // Signature: matrix(s, t, u, qt, qu, mDSqr, mQSqr)
    double c1 = 2 * (2 * p2 + omega);
    double c2 = 4 * p2 * Sin1Q * Sin2Q * std::cos(Phi1Q - Phi2Q);
    double matrixFactor = Traits::matrix(s, t, omega, qt, mDSqr, mQSqr, c1, c2);

    // Get distribution functions for p2 and p4
    double cosTheta2 = p2z / p2;
    double phi2 = std::atan2(p2y, p2x);
    double cosTheta4 = p4z / p4;
    double phi4 = std::atan2(p4y, p4x);

    double f2g = Traits::dist2(p2, cosTheta2, phi2);
    double f4g = Traits::dist4(p4, cosTheta4, phi4);

    // Statistical factor from ProcessTraits
    // We pass (0, f2, 0, f4) since f1 and f3 are not used in C(q_perp)
    double statFactor = Traits::stat(0.0, f2g, 0.0, f4g);

    // Jacobian from delta function:
    double deltaJacobian = 1.0 / (16.0 * p2 * q);

    // Full Jacobian
    double Jacobian = dOmega_dx * dp2_dx * dPhi2Q * deltaJacobian;

    double result = Jacobian * matrixFactor * statFactor;

    if (!std::isfinite(result)) {
        return 0.0;
    }

    return result;
}

} // namespace CQperpIntegral

namespace IntegrateCQperp {

std::vector<GSLVEGAS> vegasIntegrators;

/**
 * @brief Compute C(q_perp, Phi1Q) for a range of q_perp and Phi1Q values.
 *
 * @tparam ProcessTag Tag type identifying the process.
 * @param OutputFile Output file path.
 * @param p1 Momentum of particle 1 (default: 1.0).
 * @param cosTheta1 Angle of particle 1 (default: 1.0, i.e., along z-axis).
 */
template<typename ProcessTag>
void Compute(const std::string& OutputFile, double p1 = 1.0,
             double cosTheta1 = 1.0, double phi1 = 0.0) {
    auto Integrand = CQperpIntegral::CQperpIntegrand<ProcessTag>;

    size_t NqPerp = 64;
    size_t NPhi1Q = 32;
    double qPerpMin = 1e-4;
    double qPerpMax = 16.0;

    std::vector<double> qPerpValues(NqPerp);
    std::vector<double> Phi1QValues(NPhi1Q);
    std::vector<double> Results(NqPerp * NPhi1Q), Errors(NqPerp * NPhi1Q);

    for (size_t i = 0; i < NqPerp; i++) {
        qPerpValues[i] = qPerpMin * std::exp(std::log(qPerpMax / qPerpMin) * i /
                                             (NqPerp - 1));
    }

    for (size_t j = 0; j < NPhi1Q; j++) {
        Phi1QValues[j] = M_PI_2 * j / NPhi1Q;
    }

    #pragma omp parallel for collapse(2)

    for (size_t i = 0; i < NqPerp; i++) {
        for (size_t j = 0; j < NPhi1Q; j++) {
            CQperpIntegral::GSLARGS args;
            args.qPerp = qPerpValues[i];
            args.Phi1Q = Phi1QValues[j];
            args.p1 = p1;
            args.cosTheta1 = cosTheta1;
            args.phi1 = phi1;
            args.fct = Integrand;

            size_t tID = omp_get_thread_num();
            double result, error;
            vegasIntegrators[tID].integrate(args, 1e5, &result, &error);

            size_t idx = j + i * NPhi1Q;
            Results[idx] = result;
            Errors[idx] = error;
        }
    }

    // OUTPUT RESULTS //
    auto file = fmt::output_file(OutputFile);
    file.print("# q_perp  Phi1Q  C(q_perp,Phi1Q)  Error\n");

    for (size_t i = 0; i < NqPerp; i++) {
        for (size_t j = 0; j < NPhi1Q; j++) {
            size_t idx = j + i * NPhi1Q;
            file.print("{} {} {} {}\n", qPerpValues[i], Phi1QValues[j],
                       Results[idx], Errors[idx]);
        }

        file.print("\n"); // blank line for gnuplot
    }

    fmt::println(stderr, "Results written to {}", OutputFile);
}

void Setup() {
    vegasIntegrators.reserve(omp_get_max_threads());

    for (int i = 0; i < omp_get_max_threads(); i++) {
        vegasIntegrators.push_back(GSLVEGAS(CQperpIntegral::dimensions));
    }
}

} // namespace IntegrateCQperp

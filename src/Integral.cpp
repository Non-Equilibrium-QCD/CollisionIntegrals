#pragma once

#include <cmath>
#include "constants.cpp"
#include "Process.hpp"
#include <gsl/gsl_monte_vegas.h>
#include <vector>
#include "GSLVEGAS.cpp"
#include <omp.h>
#include <fmt/core.h>


namespace CollisionIntegralQCD {


/**
 * @brief Computes the collision integral for 2->2 QCD processes using Monte Carlo integration.
 * @tparam ProcessTag Tag type identifying the process (e.g., gg_to_gg, qg_to_qg).
 *         ProcessTraits<ProcessTag> must define:
 *         - dist1, dist2, dist3, dist4: distribution functions
 *         - stat: statistical factor function
 *         - matrix: matrix element squared function
 *         - nuA: degeneracy factor for particle A
 * @param x Array of random variables for sampling the integration variables.
 * @param args Struct containing additional parameters for the integration.
 * @return The value of the collision integral for the given parameters.
 */
template<typename ProcessTag>
inline double CollisionIntegral(double *x, const GSLARGS& args) {
    using Traits = ProcessTraits<ProcessTag>;

    // GET p1 ANGULAR COORDINATES //
    double p1 = args.p1;
    double cosTheta1 = args.cosTheta1;
    double Sin1 = sqrt(1.0 - cosTheta1 * cosTheta1);
    double phi1 = args.phi1;

    double p1x = p1 * Sin1 * cos(phi1);
    double p1y = p1 * Sin1 * sin(phi1);
    double p1z = p1 * cosTheta1;

    // SET BASIS VECTORS WRT p1 //
    double ep1[3]; double ep2[3]; double ep3[3];

    ep1[0] = Sin1 * cos(phi1);
    ep1[1] = Sin1 * sin(phi1);
    ep1[2] = cosTheta1;

    ep2[0] = cosTheta1 * cos(phi1);
    ep2[1] = cosTheta1 * sin(phi1);
    ep2[2] = -Sin1;

    ep3[0] = -sin(phi1);
    ep3[1] = cos(phi1);
    ep3[2] = 0.0;

    // SAMPLE p2 //
    double p2 = x[0] / (1.0 - x[0]) ; // maps [0,1] to [0,inf)
    // x[0] = p2 / (1.0 + p2);
    double Jacobian = (1.0 + p2) * (1.0 + p2);

    // SAMPLE q //
    double qMin = 0.0; double qMax = (p1 + p2);
    double q = qMin + (qMax - qMin) * x[1];
    Jacobian *= (qMax - qMin);

    // SAMPLE w //
    double wMin = std::max(-q, q - 2.0 * p2);
    double wMax = std::min(+q, 2.0 * p1 - q);
    double w = wMin + (wMax - wMin) * x[2];
    Jacobian *= (wMax - wMin);

    // SAMPLE Phi1q //
    double Phi1Q = 2.0 * M_PI * x[3];
    Jacobian *= (2.0 * M_PI);

    // SAMPLE Phi2q //
    double Phi2Q = 2.0 * M_PI * x[4];
    Jacobian *= (2.0 * M_PI);

    // GET Cos1Q AND Cos2Q FROM ENERGY CONSERVATION CONSTRAINTS //
    // double Cos1Q = -((p1 - w) * (p1 - w) - p1 * p1 - q * q) / (2.0 * p1 * q);
    double Cos1Q = w / q - (w * w - q * q) / (2.0 * p1 * q);
    double Sin1Q = sqrt(1.0 - Cos1Q * Cos1Q);
    // double Cos2Q = +((p2 + w) * (p2 + w) - p2 * p2 - q * q) / (2.0 * p2 * q);
    double Cos2Q = w / q + (w * w - q * q) / (2.0 * p2 * q);
    double Sin2Q = sqrt(1.0 - Cos2Q * Cos2Q);

    // SET q VECTOR //
    double qx = q * (Cos1Q * ep1[0] + Sin1Q * cos(Phi1Q) * ep2[0] + Sin1Q * sin(
                         Phi1Q) * ep3[0]);
    double qy = q * (Cos1Q * ep1[1] + Sin1Q * cos(Phi1Q) * ep2[1] + Sin1Q * sin(
                         Phi1Q) * ep3[1]);
    double qz = q * (Cos1Q * ep1[2] + Sin1Q * cos(Phi1Q) * ep2[2] + Sin1Q * sin(
                         Phi1Q) * ep3[2]);

    // DETERMINE ANGLES //
    double CosQ = qz / q; double SinQ = sqrt(1.0 - CosQ * CosQ);
    double PhiQ = atan2(qy, qx);

    // SET BASIS VECTORS WRT q //
    double eq1[3]; double eq2[3]; double eq3[3];

    eq1[0] = SinQ * cos(PhiQ);
    eq1[1] = SinQ * sin(PhiQ);
    eq1[2] = CosQ;

    eq2[0] = CosQ * cos(PhiQ);
    eq2[1] = CosQ * sin(PhiQ);
    eq2[2] = -SinQ;

    eq3[0] = -sin(PhiQ);
    eq3[1] = cos(PhiQ);
    eq3[2] = 0.0;

    // SET p2 VECTOR //
    double p2x = p2 * (Cos2Q * eq1[0] + Sin2Q * cos(Phi2Q) * eq2[0] + Sin2Q * sin(
                           Phi2Q) * eq3[0]);
    double p2y = p2 * (Cos2Q * eq1[1] + Sin2Q * cos(Phi2Q) * eq2[1] + Sin2Q * sin(
                           Phi2Q) * eq3[1]);
    double p2z = p2 * (Cos2Q * eq1[2] + Sin2Q * cos(Phi2Q) * eq2[2] + Sin2Q * sin(
                           Phi2Q) * eq3[2]);

    // SET p3 VECTOR //
    double p3x = p1x - qx;
    double p3y = p1y - qy;
    double p3z = p1z - qz;

    double p3 = sqrt(p3x * p3x + p3y * p3y + p3z * p3z);

    // SET p4 VECTOR //
    double p4x = p2x + qx;
    double p4y = p2y + qy;
    double p4z = p2z + qz;

    double p4 = sqrt(p4x * p4x + p4y * p4y + p4z * p4z);

    // GET MANDELSTAM VARIABLES //
    double s = +2.0 * (p1 * p2 - (p1x * p2x + p1y * p2y + p1z * p2z));
    double t = w * w - q * q;
    double u = -s - t;

    // SET qt AND qu //
    double q13 = q;
    double q23 = std::sqrt((p2x - p3x) * (p2x - p3x) + (p2y - p3y) * (p2y - p3y) +
                           (p2z - p3z) * (p2z - p3z));

    double cosTheta2 = p2z / p2;
    double phi2 = atan2(p2y, p2x);
    double cosTheta3 = p3z / p3;
    double phi3 = atan2(p3y, p3x);
    double cosTheta4 = p4z / p4;
    double phi4 = atan2(p4y, p4x);
    double f1 = Traits::dist1(p1, cosTheta1, phi1);
    double f2 = Traits::dist2(p2, cosTheta2, phi2);
    double f3 = Traits::dist3(p3, cosTheta3, phi3);
    double f4 = Traits::dist4(p4, cosTheta4, phi4);

    double s1 = Traits::stat(f1, f2, f3, f4);
    double M1 = Traits::matrix(s, t, u, q13, q23, mDSqr, mQSqr);
    double s2 = Traits::stat(f1, f2, f4, f3);
    double M2 = Traits::matrix(s, u, t, q23, q13, mDSqr, mQSqr);

    double j  = Jacobian / (16.0 * p1 * p1);
    constexpr double c  = 1.0 / power_recursive<double, 6>(2.0 * M_PI)
                          / (2.0 * Traits::nuA);

    if (std::abs(u) < 1e-12 || std::abs(t) < 1e-12) {
        return 0.0;
    }

    if (!std::isfinite(s1) || !std::isfinite(M1) || !std::isfinite(j)
            || !std::isfinite(c)) {
        fmt::println(stderr, "NaN encountered in integrand:"
                             " s1={}, M1={}, j={}, c={}",
                     s1, M1, j, c);
        return 0.0;
    }

    // return j * c * s1 * M1;
    return 0.5 * j * c * (s1 * M1 + s2 * M2);
}

}

namespace IntegrateQCD {

std::vector<GSLVEGAS> vegasIntegrators;

template<typename ProcessTag>
void Compute(double (*Integrand)(double *, const GSLARGS &)) {
    using Traits = ProcessTraits<ProcessTag>;

    size_t Np = 128;
    size_t Ncos = 16;
    double pMin = 1e-2;
    double pMax = 16.0;
    #pragma omp parallel for ordered

    for (size_t j = 0; j < Ncos; j++) {
        std::vector<double> Results(Np, 0.0), Errors(Np, 0.0), pValues(Np, 0.0);
        GSLARGS args;
        // args.cosTheta1 = 0.0;
        args.cosTheta1 = -1.0 + 2.0 * j / (Ncos - 1);
        args.phi1 = 0.0;
        args.fct = Integrand;
        size_t tID = omp_get_thread_num();

        for (size_t i = 0; i < Np; i++) {
            args.p1 = pMin * std::exp(std::log(pMax / pMin) * i / (Np - 1));
            // args.p1 = pMin + (pMax - pMin) * i / (Np - 1);

            double result, error;
            vegasIntegrators[tID].integrate(args, 10000, &result, &error);

            Results[i] = result;
            Errors[i] = error;
            pValues[i] = args.p1;
        }

        #pragma omp critical
        {
            fmt::println(stderr, "Thread {} completed ",
                         tID);
        }

        #pragma omp ordered
        {
            for (size_t i = 0; i < Np; i++) {
                // Use dist1 from traits for output (particle 1's distribution)
                fmt::println("{} {} {} {} {}", pValues[i], args.cosTheta1,
                             Traits::dist1(pValues[i], args.cosTheta1, 0.0),
                             Results[i], Errors[i]);
            }

            fmt::println("");
        }
    }

}


void Setup() {
    // Disable Cuba's internal parallelization (forking) since we use OpenMP

    vegasIntegrators.reserve(omp_get_max_threads());

    for (int i = 0; i < omp_get_max_threads(); i++) {
        vegasIntegrators.push_back(GSLVEGAS(5));
    }
}

}


// int main (int argc, const char *argv[]) {
//     IntegrateQCD::Setup();
//     IntegrateQCD::Compute();
//     return 0;
// }

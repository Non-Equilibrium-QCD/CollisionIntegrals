#pragma once

#include <cmath>
#include "constants.cpp"
#include <gsl/gsl_monte_vegas.h>
#include <vector>
#include "GSLVEGAS.cpp"
#include <omp.h>
#include <fmt/core.h>

namespace Distribution {
inline double fg(double p, double cosTheta, double phi) {
    return 1.0 / (std::expm1(p));
}

inline double fq(double p, double cosTheta, double phi) {
    return 1.0 / (std::exp(p) + 1.0);
}
}


namespace MatrixElementsQCD {

// SCREENING PARAMETER //
static const double xi0g = std::exp(5.0 / 6.0) / (2.0 * M_SQRT2);
static const double xi0q = std::exp(1.0) / M_SQRT2;

// gg -> gg //
inline double gg_gg(const double s, const double t,
                    const double u, const double qt, const double qu, const double mDSqr,
                    const double mQSqr) {

    double tBar = t * (qt * qt + xi0g * xi0g * mDSqr) / (qt * qt);
    double uBar = u * (qu * qu + xi0g * xi0g * mDSqr) / (qu * qu);

    return 4.0 * g * g * g * g * dA * CA * CA * ();

}
}

namespace StatisticalQCD {

inline double Stat(double f2, double f4) {
    return f2 * (1.0 + f4);
}
}


namespace CollisionIntegralQCD {


/**
 * @brief Computes the collision integral for 2->2 QCD processes using Monte Carlo integration.
 * @tparam StatFunc Function to compute the statistical factor.
    *         @code
    *         double StatFunc(double f1, double f2, double f3, double f4);
    *         @endcode
 * @tparam MatrixElemFunc Function to compute the squared matrix element.
    *         @code
    *         double MatrixElemFunc(double s, double t, double u, double qt, double qu, double mDSqr, double mQSqr);
    *         @endcode
 * @param x Array of random variables for sampling the integration variables.
 * @param args Struct containing additional parameters for the integration.
 * @return The value of the collision integral for the given parameters.
*/
template<auto StatFunc, auto MatrixElemFunc>
inline double CollisionIntegral(double *x, const GSLARGS& args) {

    // GET p1 ANGULAR COORDINATES //
    double qT = args.p1;
    double cosTheta1 = args.cosTheta1; double Sin1 = sqrt(1.0 - cosTheta1 * cosTheta1);
    double phi1 = args.phi1;

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

    // SAMPLE w from -oo to oo //
    double w = std::tan(M_PI * (x[2] - 0.5));
    double dw_dy = ;

    double q = std::sqrt(qT * qT + w * w);

    // SAMPLE p2 //
    double p2Min = 0.5 * (q - w);
    double p2Max = 100.0;
    double p2 = p2Min + (p2Max - p2Min) * x[0];


    // SAMPLE w //
    double wMin = std::max(-q, q - 2.0 * p2);
    double wMax = std::min(+q, 2.0 * p1 - q);
    double w = wMin + (wMax - wMin) * x[2];

    // SAMPLE Phi1q //
    double Phi1Q = 2.0 * M_PI * x[3];

    // SAMPLE Phi2q //
    double Phi2Q = 2.0 * M_PI * x[4];

    // GET Cos1Q AND Cos2Q FROM ENERGY CONSERVATION CONSTRAINTS //
    double Cos1Q = -((p1 - w) * (p1 - w) - p1 * p1 - q * q) / (2.0 * p1 * q);
    double Sin1Q = sqrt(1.0 - Cos1Q * Cos1Q);
    double Cos2Q = +((p2 + w) * (p2 + w) - p2 * p2 - q * q) / (2.0 * p2 * q);
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

    // SET p4 VECTOR //
    double p4x = p2x + qx;
    double p4y = p2y + qy;
    double p4z = p2z + qz;

    double p4 = sqrt(p4x * p4x + p4y * p4y + p4z * p4z);

    // GET MANDELSTAM VARIABLES //
    double t = w * w - q * q;

    // SET qt AND qu //
    double q13 = q;

    // SET JACOBIAN //
    double Jacobian = (2.0 * M_PI) * (2.0 * M_PI) * (p2Max - p2Min) *
                      (qMax - qMin) * (wMax - wMin) * (q * q) * (p2 * p2) * (p3 / (p1 * q)) * (p4 /
                          (p2 * q));

    double cosTheta2 = p2z / p2;
    double phi2 = atan2(p2y, p2x);
    double cosTheta4 = p4z / p4;
    double phi4 = atan2(p4y, p4x);

    double f2g = Distribution::fg(p2, cosTheta2, phi2);
    double f4g = Distribution::fg(p4, cosTheta4, phi4);

    double f2q = Distribution::fg(p2, cosTheta2, phi2);
    double f4q = Distribution::fg(p4, cosTheta4, phi4);

    double s1 = CA * f2g * (1.0 + f4g)
                + CF * f2q * (1.0 - f4q); 
    double M1 = MatrixElemFunc(s, t, u, q13, q23, mDSqr,
                mQSqr);

    double j  = Jacobian / (16.0 * p1 * p1 * p2 * p2 * q13 * q13);
    double c  = 1.0 / std::pow(2.0 * M_PI, 9);

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

    return j * c * s1 * M1;
}

}

namespace IntegrateQCD {

std::vector<GSLVEGAS> vegasIntegrators;

void Compute(double (*Integrand)(double *, const GSLARGS&)) {

    size_t Np = 64;
    size_t Ncos = 16;
    double pMin = 1e-1;
    double pMax = 4.0;
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
            // args.p1 = pMin * std::exp(std::log(pMax / pMin) * i / (Np - 1));
            args.p1 = pMin + (pMax - pMin) * i / (Np - 1);

            double result, error;
            vegasIntegrators[tID].integrate(args, 10000, &result, &error);

            Results[i] = result;
            Errors[i] = error;
            pValues[i] = args.p1;
        }

        #pragma omp critical
        {
            fmt::println(stderr, "Thread {:2} completed ",
                         fmt::styled(tID, fmt::emphasis::bold | fmt::color::green));
        }

        #pragma omp ordered
        {
            for (size_t i = 0; i < Np; i++) {
                fmt::println("{} {} {} {} {}", pValues[i], args.cosTheta1,
                             Distribution::fg(pValues[i], args.cosTheta1, 0.0),
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

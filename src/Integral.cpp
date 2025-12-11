#pragma once

#include <cmath>
#include "constants.cpp"
#include <gsl/gsl_monte_vegas.h>
#include <vector>
// #include "GSLVEGAS.cpp"
#include "CUBAVEGAS.cpp"
#include <omp.h>
#include <fmt/core.h>


struct angles {
    double cosTheta, sinTheta, phi;
};

class vec3 {
public:
    double x, y, z;
    double cosTheta, sinTheta, phi;
    vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {
        double r = std::sqrt(x * x + y * y + z * z);
        cosTheta = z / r;
        sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
        phi = std::atan2(y, x);
    }

    vec3(const angles &Ref, double q, double cosTheta1q, double sinTheta1q,
         double phiq) {

        ///////////////////////////
        // SET COORDINATE SYSTEM //
        ///////////////////////////

        if (std::abs(Ref.sinTheta) < 1e-12) {
            // Ref along z //
            x = q * sinTheta1q * std::cos(phiq);
            y = q * sinTheta1q * std::sin(phiq);
            z = q * cosTheta1q;

            // if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
            //     fmt::print(stderr, "NaN encountered in vec3 constructor\n");
            // }

            cosTheta = cosTheta1q;
            sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
            phi = phiq;
            return;

        }

        // e1=eP1 //
        double e1x = Ref.sinTheta * std::cos(Ref.phi);
        double e1y = Ref.sinTheta * std::sin(Ref.phi);
        double e1z = Ref.cosTheta;

        // e2=eZ-CosThetaQ eQ //
        double e2x = -Ref.cosTheta * e1x;
        double e2y = -Ref.cosTheta * e1y;
        double e2z = 1.0 - Ref.cosTheta * e1z;

        double e2Abs = std::sqrt(e2x * e2x + e2y * e2y + e2z * e2z);

        e2x /= e2Abs; e2y /= e2Abs; e2z /= e2Abs;

        // e3=e1 x e2 //
        double e3x = e1y * e2z - e2y * e1z;
        double e3y = e1z * e2x - e2z * e1x;
        double e3z = e1x * e2y - e2x * e1y;

        double e3Abs = std::sqrt(e3x * e3x + e3y * e3y + e3z * e3z);

        e3x /= e3Abs; e3y /= e3Abs; e3z /= e3Abs;

        x = q * (cosTheta1q * e1x + sinTheta1q * std::cos(phiq) * e2x
                 + sinTheta1q * std::sin(phiq) * e3x);
        y = q * (cosTheta1q * e1y + sinTheta1q * std::cos(phiq) * e2y
                 + sinTheta1q * std::sin(phiq) * e3y);
        z = q * (cosTheta1q * e1z + sinTheta1q * std::cos(phiq) * e2z
                 + sinTheta1q * std::sin(phiq) * e3z);

        // if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
        //     fmt::print(stderr, "NaN encountered in vec3 constructor\n");
        // }

        cosTheta = z / q;
        sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
        phi = std::atan2(y, x);
    }
};

namespace Distribution {
inline double fg(double p, double cosTheta, double phi) {
    // return std::exp(- (10.0 - p) * (10.0 - p));
    return std::exp(- p * p) * std::exp(- cosTheta * cosTheta);
}
}


namespace MatrixElementsQCD {

// SCREENING PARAMETER //
static const double xi0g = std::exp(5.0 / 6.0) / (2.0 * M_SQRT2);
static const double xi0q = std::exp(1.0) / M_SQRT2;

// gg -> gg //
inline double gg_gg_MatrixElementSqr(const double s, const double t,
                                     const double u, const double qt, const double qu, const double mDSqr,
                                     const double mQSqr) {

    double tBar = t * (qt * qt + xi0g * xi0g * mDSqr) / (qt * qt);
    double uBar = u * (qu * qu + xi0g * xi0g * mDSqr) / (qu * qu);

    double uchannel = (s - t) * (s - t) / (uBar * uBar);
    double tchannel = (s - u) * (s - u) / (tBar * tBar);

    if (abs(uchannel) > 1e10) {
        fmt::println(stderr, "s={}, t={}, u={}, tBar={}, uBar={}, uchannel={}",
                     s, t, u, tBar, uBar, uchannel);
    }

    if (abs(tchannel) > 1e10) {
        fmt::println(stderr, "s={}, t={}, u={}, tBar={}, uBar={}, tchannel={}",
                     s, t, u, tBar, uBar, tchannel);
    }

    return 4.0 * g * g * g * g * dA * CA * CA * (9.0
            + (s - u) * (s - u) / (tBar * tBar)
            + (s - t) * (s - t) / (uBar * uBar)
            + (u - t) * (u - t) / (s * s));

}
}

namespace StatisticalQCD {

inline double gg_gg(double f1, double f2, double f3, double f4) {
    return f1 * f2 * (1.0 + f3) * (1.0 + f4) - f3 * f4 * (1.0 + f1) * (1.0 + f2);
}
}


namespace CollisionIntegralQCD {
inline double gg_gg_CollisionIntegral(double *x, const GSLARGS& args) {

    double p1 = args.p1;
    double cosTheta1 = args.cosTheta1;
    double sinTheta1 = std::sqrt(1.0 - cosTheta1 * cosTheta1);
    double phi1 = args.phi1;

    // p2 Integral //
    double p2 = std::tan(M_PI_2 * x[0]); // Map [0,1] to [0, infinity)
    double Jacobian = M_PI_2 * (1.0 + p2 * p2) * p2 * p2;
    double phi2 = 2.0 * M_PI * x[1]; // Map [0,1] to [0,2pi]
    Jacobian *= 2.0 * M_PI;

    // q Integral //
    double qMin = 0.0;
    double qMax = p1 + p2;
    double qt = qMin + (qMax - qMin) * x[2]; // Map [0,1] to [qMin, qMax]
    Jacobian *= (qMax - qMin) * qt * qt;
    double phiq = 2.0 * M_PI * x[3]; // Map [0,1] to [0,2pi]
    Jacobian *= 2.0 * M_PI;

    // max(-q, -2*p2 + q) < w, w < min(2*p1 - q, q)
    double wMin = std::max(-2.0 * p2 + qt, -qt);
    double wMax = std::min(2.0 * p1 - qt, qt);

    if (wMax <= wMin) { return 0.0; }

    double w = wMin + (wMax - wMin) * x[4];
    Jacobian *= (wMax - wMin);
    double t = w * w - qt * qt;

    double cosTheta1q = w / qt - t / (2.0 * p1 * qt);
    double cosTheta2q = w / qt + t / (2.0 * p2 * qt);
    double sinTheta1q = std::sqrt(1.0 - cosTheta1q * cosTheta1q);
    double sinTheta2q = std::sqrt(1.0 - cosTheta2q * cosTheta2q);


    double cosTheta12 = cosTheta1q * cosTheta2q
                        + sinTheta1q * sinTheta2q * std::cos(phi1 - phi2);
    double s = 2.0 * p1 * p2 * (1.0 - cosTheta12);
    double u = -s - t;

    double p1x = p1 * sinTheta1 * std::cos(phi1);
    double p1y = p1 * sinTheta1 * std::sin(phi1);
    double p1z = p1 * cosTheta1;

    angles Ref1 = {cosTheta1, sinTheta1, phi1};
    vec3 qVec(Ref1, qt, cosTheta1q, sinTheta1q, phiq);
    double qx = qVec.x;
    double qy = qVec.y;
    double qz = qVec.z;


    angles Ref2 = {qVec.cosTheta, qVec.sinTheta, qVec.phi};
    vec3 p2Vec(Ref2, p2, cosTheta2q, sinTheta2q, phi2);
    double cosTheta2 = p2Vec.cosTheta;
    double p2x = p2Vec.x;
    double p2y = p2Vec.y;
    double p2z = p2Vec.z;

    double p3x = p1x - qx;
    double p3y = p1y - qy;
    double p3z = p1z - qz;

    double p4x = p2x + qx;
    double p4y = p2y + qy;
    double p4z = p2z + qz;

    // (p2 - p3)^2 = (p2 - p1 + q)^2
    double qu = std::sqrt(
                    (p2x - p3x) * (p2x - p3x)
                    + (p2y - p3y) * (p2y - p3y)
                    + (p2z - p3z) * (p2z - p3z)
                );


    double p3 = std::sqrt(p3x * p3x + p3y * p3y + p3z * p3z);
    double p4 = std::sqrt(p4x * p4x + p4y * p4y + p4z * p4z);
    double cosTheta3 = p3z / p3;
    double phi3 = std::atan2(p3y, p3x);
    double cosTheta4 = p4z / p4;
    double phi4 = std::atan2(p4y, p4x);

    double f1 = Distribution::fg(p1, cosTheta1, phi1);
    double f2 = Distribution::fg(p2, cosTheta2, phi2);
    double f3 = Distribution::fg(p3, cosTheta3, phi3);
    double f4 = Distribution::fg(p4, cosTheta4, phi4);

    double s1 = StatisticalQCD::gg_gg(f1, f2, f3, f4);
    double M1 = MatrixElementsQCD::gg_gg_MatrixElementSqr(s, t, u, qt, qu, mDSqr,
                mQSqr);
    // double s2 = StatisticalQCD::gg_gg(f1, f2, f4, f3);
    // double M2 = MatrixElementsQCD::gg_gg_MatrixElementSqr(s, u, t, qu, qt, mDSqr,
    //             mQSqr);

    double j  = Jacobian / (16.0 * p1 * p1 * p2 * p2 * qt * qt);
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

    if (abs(s1 * M1 * j * c) > 1e2) {
        fmt::println(stderr, "{} \n", s1 * M1 * j * c);
    }

    // fmt::println(stderr, "{}", qu);
    return j * c * s1 * M1;
    // return j * c * 0.5 * (s1 * M1 + s2 * M2);
}

}

namespace IntegrateQCD {

std::vector<GSLVEGAS> vegasIntegrators;

void Compute() {

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
        args.fct = CollisionIntegralQCD::gg_gg_CollisionIntegral;
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
            fmt::println(stderr, "Thread {} completed ",
                         tID);
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
    int ncores = 0;
    int pcores = 0;
    cubacores(&ncores, &pcores);

    vegasIntegrators.reserve(omp_get_max_threads());

    for (int i = 0; i < omp_get_max_threads(); i++) {
        vegasIntegrators.push_back(GSLVEGAS(5));
    }
}

}

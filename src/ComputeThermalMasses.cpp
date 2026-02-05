#pragma once
#pragma once
#include <array>
#include <cmath>
#include "constants.cpp"
#include <gsl/gsl_monte_vegas.h>
#include <vector>
#include "GSLVEGAS.cpp"
#include "GAUSSQUADRATURE.cpp"
#include <omp.h>
#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/color.h>

typedef QUAD<3, 52> GaussQuad3D;

namespace ThermalMasses {
struct GSLARGS {
    double (*fct)(double *, const GSLARGS &);
};

std::vector<GSLVEGAS> vegasIntegrators;
constexpr size_t dimensions = 3;

template<double (* Dist_g)(double, double, double)>
inline double CollisionIntegral(double *x, const GSLARGS& args) {

    // GET p1 ANGULAR COORDINATES //
    // SAMPLE p2 //
    double p1 = x[0] / (1.0 - x[0]) ; // maps [0,1] to [0,inf)

    if (x[0] == 1.0) {
        return 0.0;
    }

    // x[0] = p2 / (1.0 + p2);
    double cosTheta1 = -1.0 + 2.0 * x[1];
    double phi1 = 2.0 * M_PI * x[2];
    double Jacobian = (1.0 + p1) * (1.0 + p1) * (2.0 * M_PI) * 2.0;

    double f1 = Dist_g(p1, cosTheta1, phi1);

    double j  = Jacobian / (2.0 * p1);

    double Measure = p1 * p1 / std::pow(2.0 * M_PI, 3);

    if (!std::isfinite(f1) || !std::isfinite(j) ) {
        fmt::println(stderr, "NaN encountered in integrand:"
                             " f1={}, j={}",
                     f1, j);
        return 0.0;
    }

    return j * Measure * ( f1 );
}

template<double (* Dist_g)(double, double, double), double (* Dist_q)(double, double, double)>
std::array<double, 2> Compute() {

    GSLARGS args;
    args.fct = CollisionIntegral<Dist_g>;;
    size_t tID = omp_get_thread_num();
    double r1, e1, r2, e2;
    vegasIntegrators[tID].integrate(args, static_cast<int>(1e6), &r1, &e1);
    args.fct = CollisionIntegral<Dist_q>;
    vegasIntegrators[tID].integrate(args, static_cast<int>(1e7), &r2, &e2);

    return {r1, r2};

}

double mDfct(double Ig, double Iq) {
    return 4.0 * g * g * ( 2.0 * Nc * Ig + 2.0 * Nf * Iq);
}

double mQfct(double Ig, double Iq) {
    return g * g * CF * ( 2.0 * Ig + 2.0 * Iq);
}
void Setup() {
    vegasIntegrators.reserve(omp_get_max_threads());

    for (int i = 0; i < omp_get_max_threads(); i++) {
        vegasIntegrators.push_back(GSLVEGAS(ThermalMasses::dimensions));
    }
}
};

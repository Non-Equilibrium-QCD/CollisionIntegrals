#pragma once

#include <cmath>
#include <gsl/gsl_monte_vegas.h>
#include <vector>
#include "GSLVEGAS.cpp"
#include <omp.h>
#include <fmt/core.h>


inline double fg(double p) {
    return std::exp(- p * p);
}


inline double Integral(double *x, const GSLARGS& args) {

    // // p2 Integral //
    double p2 = std::tan(M_PI_2 * x[0]); // Map [0,1] to [0, infinity)
    double Jacobian = M_PI_2 * (1.0 + p2 * p2);
    return fg(p2) * Jacobian;
}

namespace IntegrateQCD {

std::vector<GSLVEGAS> vegasIntegrators;

void Compute() {

    size_t tID = omp_get_thread_num();
    GSLARGS args;
    args.fct = Integral;

    double result, error;
    vegasIntegrators[tID].integrate(args, 1000, &result, &error);

    fmt::println("Integral: {} Error: {}", result, error);
    fmt::println("Exact: {}", 0.5 * std::sqrt(M_PI));

}


void Setup() {
    vegasIntegrators.reserve(omp_get_max_threads());

    for (int i = 0; i < omp_get_max_threads(); i++) {
        vegasIntegrators.push_back(GSLVEGAS(1));
    }
}

}


int main (int argc, char *argv[]) {
    IntegrateQCD::Setup();
    IntegrateQCD::Compute();
    return 0;
}

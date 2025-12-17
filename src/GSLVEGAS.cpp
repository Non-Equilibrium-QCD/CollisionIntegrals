#pragma once

#include <cmath>
#include <cstddef>
#include <gsl/gsl_monte_vegas.h>
#include <vector>

struct GSLARGS {
    double p1, cosTheta1, phi1;
    double (*fct)(double *, const GSLARGS &);
};

namespace private_GSLVEGAS {

double Integrand(double *x, size_t dim, void *params) {
    (void)dim;  // Unused parameter
    GSLARGS args = *static_cast<GSLARGS *>(params);
    return args.fct(x, args);
}

}

class GSLVEGAS {
public:
    explicit GSLVEGAS(size_t dim): dimension(dim) {
        // vegasState = gsl_monte_vegas_alloc(dimension);
        rng = gsl_rng_alloc(gsl_rng_default);
    }

    // Delete copy (raw pointers can't be safely copied)
    GSLVEGAS(const GSLVEGAS &) = delete;
    GSLVEGAS &operator=(const GSLVEGAS &) = delete;

    // Move constructor
    GSLVEGAS(GSLVEGAS&& other) noexcept
        : dimension(other.dimension), /* vegasState(other.vegasState),*/ rng(other.rng) {
        // other.vegasState = nullptr;
        other.rng = nullptr;
    }

    // Move assignment
    GSLVEGAS &operator=(GSLVEGAS&& other) noexcept {
        if (this != &other) {
            // if (vegasState) { gsl_monte_vegas_free(vegasState); }

            if (rng) { gsl_rng_free(rng); }

            dimension = other.dimension;
            // vegasState = other.vegasState;
            rng = other.rng;
            // other.vegasState = nullptr;
            other.rng = nullptr;
        }

        return *this;
    }

    void integrate(const GSLARGS& args,
                   size_t calls, double *result, double *error) {
        std::vector<double> xl(dimension, 0.0);
        std::vector<double> xu(dimension, 1.0);
        gsl_monte_function f;
        f.f = &private_GSLVEGAS::Integrand;
        f.dim = dimension;
        f.params = const_cast<GSLARGS *>(&args);
        auto vs = gsl_monte_vegas_alloc(dimension);
        gsl_monte_vegas_integrate(&f, xl.data(), xu.data(), dimension, calls, rng,
                                  vs, result, error);
        size_t iteration = 0;

        while (std::fabs(gsl_monte_vegas_chisq(vs) - 1.0) > 0.5 and
                iteration < 20) {
            gsl_monte_vegas_integrate(&f, xl.data(), xu.data(), dimension, calls / 5, rng,
                                      vs, result, error);
            iteration++;
        }
    }

    ~GSLVEGAS() {
        // if (vegasState) { gsl_monte_vegas_free(vegasState); }

        if (rng) { gsl_rng_free(rng); }
    }

private:
    size_t dimension;
    // gsl_monte_vegas_state *vegasState;
    gsl_rng *rng;
};

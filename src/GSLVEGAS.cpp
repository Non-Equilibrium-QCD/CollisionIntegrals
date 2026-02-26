#pragma once

#include <cmath>
#include <cstddef>
#include <gsl/gsl_monte_vegas.h>
#include <vector>

namespace private_GSLVEGAS {

template<typename ARGTYPE>
double Integrand(double *x, size_t dim, void *params) {
    (void)dim;  // Unused parameter
    ARGTYPE args = *static_cast<ARGTYPE *>(params);
    return args.fct(x, args);
}

}

class GSLVEGAS {
public:
    explicit GSLVEGAS(size_t dim): dimension(dim) {
        rng = gsl_rng_alloc(gsl_rng_default);
    }

    // Delete copy (raw pointers can't be safely copied)
    GSLVEGAS(const GSLVEGAS &) = delete;
    GSLVEGAS &operator=(const GSLVEGAS &) = delete;

    // Move constructor
    GSLVEGAS(GSLVEGAS&& other) noexcept
        : dimension(other.dimension), rng(other.rng) {
        other.rng = nullptr;
    }

    // Move assignment
    GSLVEGAS &operator=(GSLVEGAS&& other) noexcept {
        if (this != &other) {

            if (rng) { gsl_rng_free(rng); }

            dimension = other.dimension;
            rng = other.rng;
            other.rng = nullptr;
        }

        return *this;
    }

    template<typename ARGTYPE>
    void integrate(const ARGTYPE& args,
                   size_t calls, double *result, double *error) {
        std::vector<double> xl(dimension, 0.0);
        std::vector<double> xu(dimension, 1.0);
        gsl_monte_function f;
        f.f = &private_GSLVEGAS::Integrand<ARGTYPE>;
        f.dim = dimension;
        f.params = const_cast<ARGTYPE *>(&args);
        auto vs = gsl_monte_vegas_alloc(dimension);
        
        // Warm-up phase
        size_t warmup_calls = calls / 10;
        if (warmup_calls < 10000) warmup_calls = 10000;
        
        for (size_t i = 0; i < 5; ++i) {
            gsl_monte_vegas_integrate(&f, xl.data(), xu.data(), dimension, warmup_calls, rng,
                                      vs, result, error);
        }
        
        // Main integration phase
        gsl_monte_vegas_integrate(&f, xl.data(), xu.data(), dimension, calls, rng,
                                  vs, result, error);
                                  
        size_t iteration = 0;
        while (std::fabs(gsl_monte_vegas_chisq(vs) - 1.0) > 0.5 &&
                iteration < 10) {
            gsl_monte_vegas_integrate(&f, xl.data(), xu.data(), dimension, calls, rng,
                                      vs, result, error);
            iteration++;
        }
        
        gsl_monte_vegas_free(vs);
    }

    ~GSLVEGAS() {

        if (rng) { gsl_rng_free(rng); }
    }

private:
    size_t dimension;
    gsl_rng *rng;
};

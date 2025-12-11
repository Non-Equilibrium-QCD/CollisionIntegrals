#pragma once

// CUBA DEFINITIONS //
#define MINEVAL 3200
#define MAXEVAL 2560000

// VEGAS PARAMETERS //
#define NSTART 10
#define NINCREASE 500
#define NBATCH 100
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

// SUAVE PARAMETERS //
#define NNEW 8000
#define NMIN 8
#define FLATNESS 15.

/// CUHRE PARAMETERS //
#define KEY 0



// DIVONNE PARAMETERS //
#define KEY1 11
#define KEY2 11
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0


#include <cmath>
#include <cstddef>
#include <cuba.h>
#include <vector>

struct GSLARGS {
    double p1, cosTheta1, phi1;
    double (*fct)(double *, const GSLARGS &);
};

namespace private_GSLVEGAS {

static int Integrand(const int *ndim, const double x[],
                     const int *ncomp, double fval[], void *userdata) {
    GSLARGS args = *static_cast<GSLARGS *>(userdata);
    fval[0] = args.fct(const_cast<double *>(x), args);
    return 1;
}

}

class GSLVEGAS {
public:
    explicit GSLVEGAS(size_t dim): dimension(dim) {
    }

    // Delete copy (raw pointers can't be safely copied)
    GSLVEGAS(const GSLVEGAS &) = delete;
    GSLVEGAS &operator=(const GSLVEGAS &) = delete;

    // Move constructor
    GSLVEGAS(GSLVEGAS&& other) noexcept
        : dimension(other.dimension) {
    }

    // Move assignment
    GSLVEGAS &operator=(GSLVEGAS&& other) noexcept {
        if (this != &other) {
            dimension = other.dimension;
        }

        return *this;
    }

    void integrate(const GSLARGS& args,
                   size_t calls, double *result, double *error) {
        int nregions, neval, fail;
        double prob;
        double epsrel = 1e-5; double epsabs = 1e-5;
        // Cuhre(dimension, 1, private_GSLVEGAS::Integrand, const_cast<GSLARGS *>(&args),
        //       1,
        //       epsrel, epsabs, 0,
        //       MINEVAL, MAXEVAL, 0,
        //       STATEFILE, SPIN,
        //       &nregions, &neval, &fail, result, error, &prob);

        // Vegas(dimension, 1, private_GSLVEGAS::Integrand, const_cast<GSLARGS *> (&args),
        //       1,
        //       epsrel, epsabs, 0, 42,
        //       MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
        //       GRIDNO, STATEFILE, SPIN,
        //       &neval, &fail, result, error, &prob);

        Suave(dimension, 1, private_GSLVEGAS::Integrand, const_cast<GSLARGS *> (&args),
              1,
              epsrel, epsabs, 0, 42,
              MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
              STATEFILE, SPIN,
              &nregions, &neval, &fail, result, error, &prob);

        // Divonne(dimension, 1, private_GSLVEGAS::Integrand,
        //         const_cast<GSLARGS *> (&args), 1,
        //         epsrel, epsabs, 0, 42,
        //         MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
        //         BORDER, MAXCHISQ, MINDEVIATION,
        //         NGIVEN, dimension, NULL, NEXTRA, NULL,
        //         STATEFILE, SPIN,
        //         &nregions, &neval, &fail, result, error, &prob);
        
    }

    ~GSLVEGAS() {
    }

private:
    size_t dimension;
};

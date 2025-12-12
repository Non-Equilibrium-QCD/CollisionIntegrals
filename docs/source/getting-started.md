# Getting Started

## Prerequisites

- C++17 compiler (GCC 7+, Clang 5+)
- GSL (GNU Scientific Library)
- OpenMP
- fmt library

## Building

```bash
make          # build main executable
make testint  # build test executable
make docs     # generate documentation
make docs-pdf # generate PDF documentation
```

## Quick Example

We use template aliases to pass the statistical factor and matrix element to the collision integral framework.

```cpp
#include "Integral.cpp"

namespace StatisticalQCD {
inline double gg_gg(double f1, double f2, double f3, double f4) {
    return f1 * f2 * (1.0 + f3) * (1.0 + f4) - f3 * f4 * (1.0 + f1) * (1.0 + f2);
}
}


namespace MatrixElementsQCD {
inline double gg_gg_MatrixElementSqr(const double s, const double t,
                                     const double u, const double qt, const double qu, const double mDSqr,
                                     const double mQSqr) {
    double tBar = t * (qt * qt + xi0g * xi0g * mDSqr) / (qt * qt);
    double uBar = u * (qu * qu + xi0g * xi0g * mDSqr) / (qu * qu);

    return 4.0 * g * g * g * g * dA * CA * CA * (9.0
            + (s - u) * (s - u) / (tBar * tBar)
            + (s - t) * (s - t) / (uBar * uBar)
            + (u - t) * (u - t) / (s * s));

}
}

int main() {
    IntegrateQCD::Setup();
    auto gg_gg = CollisionIntegralQCD::CollisionIntegral<
        StatisticalQCD::gg_gg,
        MatrixElementsQCD::gg_gg_MatrixElementSqr>;
    IntegrateQCD::Compute(gg_gg);
    return 0;
}
```

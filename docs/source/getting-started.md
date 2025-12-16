# Getting Started

## Prerequisites

- C++17 compiler (GCC 7+, Clang 5+)
- GSL (GNU Scientific Library)
- Cuba library
- OpenMP
- fmt library

## Building

```bash
make          # build main executable
make test     # run all tests
make docs     # generate HTML documentation
make docs-pdf # generate PDF documentation
```

## Examples

Example programs are located in the `examples/` directory. Build and run them with:

```bash
make example NAME=YM       # build examples/YM.cpp
make example-run NAME=YM   # build and run
```

Running `make example` without `NAME=` shows available examples:

```
Usage: make example NAME=<name>
       make example-run NAME=<name>

Available examples:
  - YM

Example: make example NAME=YM
```

### Creating Your Own Example

Create a new `.cpp` file in `examples/`:

```cpp
// examples/MyExample.cpp
#include "../src/Integral.cpp"

int main() {
    IntegrateQCD::Setup();
    auto integrand = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>;
    IntegrateQCD::Compute<gg_to_gg>(integrand);
    return 0;
}
```

Then build with `make example NAME=MyExample`.

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

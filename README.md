# CollIntegral

QCD Collision Integral Computation Library for elastic 2 -> 2 scattering processes.

## Overview

This library computes collision integrals for QCD processes using Monte Carlo integration. It supports multiple integration backends (GSL VEGAS, Cuba SUAVE/VEGAS/Cuhre) with OpenMP parallelization.

## Prerequisites

- C++17 compiler (GCC 7+, Clang 5+)
- [GSL](https://www.gnu.org/software/gsl/) (GNU Scientific Library)
- [Cuba](http://www.feynarts.de/cuba/) library
- [fmt](https://fmt.dev/) library
- OpenMP

### Ubuntu/Debian

```bash
sudo apt install libgsl-dev libfmt-dev
```

For Cuba, download and install from: http://www.feynarts.de/cuba/

## Building

```bash
make          # build main executable
make testint  # build test executable
make run      # build and run
make debug    # build with debug flags
```

## Documentation

Documentation is generated using Doxygen with the doxygen-awesome theme.

```bash
make docs          # generate HTML docs
make docs-open     # generate and open in browser
make docs-pdf      # generate PDF (requires texlive)
make docs-pdf-open # generate PDF and open
make docs-clean    # remove generated docs
```

## Project Structure

```
src/
  main.cpp          # Entry point
  Integral.cpp      # Collision integral implementation
  GSLVEGAS.cpp      # GSL VEGAS integration wrapper
  CUBAVEGAS.cpp     # Cuba library integration wrapper
  constants.cpp     # Physical constants
  QCD/
    GLUON.cpp       # Gluon-specific calculations
docs/
  pages/            # Documentation source files
```

## Usage

```cpp
#include "Integral.cpp"

int main() {
    IntegrateQCD::Setup();
    IntegrateQCD::Compute();
    return 0;
}
```

## License

[Add license here]

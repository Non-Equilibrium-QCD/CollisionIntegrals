# Collision Integral Library

Collision Integral Computation Library for elastic 2 <-> 2 elastic scattering processes.

## Overview

This library computes collision integrals for QCD processes using Monte Carlo integration. It supports GSL VEGAS with OpenMP parallelization.

## Prerequisites

- C++17 compiler: recent GCC
- [GSL](https://www.gnu.org/software/gsl/) (GNU Scientific Library)
- [fmt](https://fmt.dev/) library (Included in `third_party/`)
- [doctest](https://github.com/doctest/doctest) library (Included in `third_party/`)
- OpenMP


## Usage

Compile and run the YM example:
```bash
make example NAME=YM       # build examples/YM.cpp
```
This will compute the collision integral for gluon-gluon scattering in `OUTPUT/QCDgg_gg.dat`.
Take a look at the docs for how to compute your own processes!

## Documentation

Documentation is generated using Doxygen with the doxygen-awesome theme.

```bash
make docs          # generate HTML docs
make docs-open     # generate and open in browser
make docs-pdf      # generate PDF (requires texlive)
make docs-pdf-open # generate PDF and open
make docs-clean    # remove generated docs
```

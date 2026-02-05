(FunctionExpansion)=
# Function Expansion
In this section, we describe `FunctionExpansion.cpp`, which provides tools to decompose the phase-space distribution into `orthogonal` basis functions.
The phase-space distribution $f(p)$ can be expanded as:
```{math}
:label: FunctionExpansion
f(p, \cos\theta, \phi) = \sum_{ijk} a_{ijk} R_i(p) Y_{jk}(\cos\theta, \phi)
```

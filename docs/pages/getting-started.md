# Getting Started {#getting-started}

## Prerequisites

- C++17 compiler (GCC 7+, Clang 5+)
- GSL (GNU Scientific Library)
- OpenMP
- fmt library

## Phase-Space Measure

The collision integral for \f$ 2 \to 2 \f$ scattering processes is given by:

\f{align}{
    C[f_1] = \frac{1}{2E_1} \int \frac{d^3p_2}{(2\pi)^3 2E_2} 
    \frac{d^3p_3}{(2\pi)^3 2E_3} \frac{d^3p_4}{(2\pi)^3 2E_4} 
    (2\pi)^4 \delta^{(4)}(p_1 + p_2 - p_3 - p_4) \, 
    |\mathcal{M}|^2 \, \mathcal{F}[f]\;,
\f}
where the statistical factor is:
\f{align}{
    \mathcal{F}[f] = f_1 f_2 (1 \pm f_3)(1 \pm f_4) - f_3 f_4 (1 \pm f_1)(1 \pm f_2)
\f}


We first introduce the momentum transfer \f$ q = p_1 - p_3 \f$ and the energy transfer \f$ \omega = E_1 - E_3 \f$.
Writing the energy delta function as:
\f{align}{
    \delta(E_1 + E_2 - E_3 - E_4) = \int d\omega \, \delta(E_1 - \omega - E_3) \, \delta(E_2 + \omega - E_4)
\f}
Using the delta functions we perform the integral over \f$ p_4 \f$ and shift the integral \f$ p_3 \to q \f$.
The positivity of energies \f$ E_3, E_4 \geq 0 \f$ implies \f$ -q \leq \omega \leq q \f$.
We obtain the reduced collision integral:
\f{align}{
    C[f_1] = & \int \frac{d^3p_2}{(2\pi)^3}
    \frac{d^3 q}{(2\pi)^3}
    \int_{-q}^{q} \frac{d\omega}{2\pi} \,
    \frac{1}{16p_1 p_2 |\vec{p}_1 - \vec{q}| |\vec{p}_2 + \vec{q}|} \\
    & \times (2\pi) \delta(p_1 - \omega - |\vec{p}_1 - \vec{q}|) \,
    (2\pi) \delta(p_2 + \omega - |\vec{p}_2 + \vec{q}|) \,
    |\mathcal{M}|^2 \, \mathcal{F}[f]\;.
\f}

We use the angles \f$\cos\theta_{1q} = (\vec{p}_1 \cdot \vec{q})/(p_1 q)\f$ and \f$\cos\theta_{2q} = (\vec{p}_2 \cdot \vec{q})/(p_2 q)\f$ to perform the integrals over the delta functions.
The Jacobians are:
\f{align}{
   \delta(p_1 - \omega - |\vec{p}_1 - \vec{q}|)
    &= \frac{|\vec{p}_1 - \vec{q}|}{p_1 q} \,
     \delta\left(\cos\theta_{1q} - \left(\frac{\omega}{q} - \frac{\omega^2 - q^2}{2 p_1 q}\right)\right)
     \Theta\left(p_1 - \frac{q + \omega}{2}\right)\;,\\
   \delta(p_2 + \omega - |\vec{p}_2 + \vec{q}|)
    &= \frac{|\vec{p}_2 + \vec{q}|}{p_2 q} \,
     \delta\left(\cos\theta_{2q} - \left(\frac{\omega}{q} + \frac{\omega^2 - q^2}{2 p_2 q}\right)\right)
     \Theta\left(p_2 - \frac{q - \omega}{2}\right)\;,
\f}
where the \f$\Theta\f$-functions together with the positivity of \f$p_1\f$ and \f$p_2\f$ ensures the conditions \f$ |\cos\theta_{1q}| \leq 1 \f$ and \f$ |\cos\theta_{2q}| \leq 1 \f$.

Combining the \f$\Theta\f$-functions, we obtain the integration limits:
\f{align}{
    q &\leq p_1 + p_2\;, \\
    \max(-q, q - 2p_2) &\leq \omega \leq \min(q, 2p_1 - q)\;.
\f}

Finally, we obtain the expression for the collision integral:
\f{align}{
    C[f_1] = & \frac{1}{{(2\pi)}^6}
    \int dp_2 \int_{0}^{p_1+p_2} dq
    d\cos\theta_{2q} d\phi_{2q}
    d\cos\theta_{1q} d\phi_{1q}
    \int_{\min(q, 2p_1 - q)}^{\max(-q, q - 2p_2)} \frac{d\omega}{2\pi} \,
    \frac{p_2^2 q^2}{16p_1^2 p_2^2 q^2} \\
    & \times (2\pi) \delta\left(\cos\theta_{1q} - \left(\frac{\omega}{q} - \frac{\omega^2 - q^2}{2 p_1 q}\right)\right) \,
    (2\pi) \delta\left(\cos\theta_{2q} - \left(\frac{\omega}{q} + \frac{\omega^2 - q^2}{2 p_2 q}\right)\right) \,
    |\mathcal{M}|^2 \, \mathcal{F}[f]\;.
\f}

We introduce the basis vectors:
\f{align}{
    \vec{e}^{\, 1}_{p_1} &= (\sin\theta_{p_1} \cos\phi_{p_1}, \sin\theta_{p_1} \sin\phi_{p_1}, \cos\theta_{p_1})\;, \\
    \vec{e}^{\, 2}_{p_1} &= (\cos\theta_{p_1} \cos\phi_{p_1}, \cos\theta_{p_1} \sin\phi_{p_1}, -\sin\theta_{p_1})\;, \\
    \vec{e}^{\,3}_{p_1} &= (-\sin\phi_{p_1}, \cos\phi_{p_1}, 0)\;.
\f}
Allowing us to express the momentum \f$ \vec{q}\f$ as:
\f{align}{
    \vec{q} = q \left[
     \cos\theta_{1q}\vec{e}^{\, 1}_{p_1}
     +\sin\theta_{1q} \left(\cos\phi_{1q}\vec{e}^{\, 2}_{p_1} + \sin\phi_{1q}\vec{e}^{\, 3}_{p_1}\right)
      \right]\;.
\f}

Similarly, we can express \f$ \vec{p}_2 \f$ in terms of the basis vectors defined by \f$ \vec{q} \f$:
\f{align}{
    \vec{e}^{\, 1}_{q} &= (\sin\theta_{q} \cos\phi_{q}, \sin\theta_{q} \sin\phi_{q}, \cos\theta_{q})\;, \\
    \vec{e}^{\, 2}_{q} &= (\cos\theta_{q} \cos\phi_{q}, \cos\theta_{q} \sin\phi_{q}, -\sin\theta_{q})\;, \\
    \vec{e}^{\,3}_{q} &= (-\sin\phi_{q}, \cos\phi_{q}, 0)\;.
\f}
Leading to the expression:
\f{align}{
    \vec{p}_2 = p_2 \left[
     \cos\theta_{2q}\vec{e}^{\, 1}_{q}
     +\sin\theta_{2q} \left(\cos\phi_{2q}\vec{e}^{\, 2}_{q} + \sin\phi_{2q}\vec{e}^{\, 3}_{q}\right)
      \right]\;.
\f}



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

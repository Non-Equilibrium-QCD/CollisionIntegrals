(SimpleIntegral)=
# Test Integral

Starting from the final form of the collision integral in {eq}`eq-final-collision-integral`, we compute the collision integral for constant matrix element $|\mathcal{M}|^2 = 1$.
We take the statistical term as a simple Gaussian for $p_2$ as $f(p_2) = e^{-p^2}$ leading to the following form of the collision integral:
```{math}
:label:
C[f_1] = \frac{1}{2\nu_1}\frac{1}{{(2\pi)}^4}\frac{1}{16p_1^2}
\int_0^\infty dp_2 \int_{0}^{p_1+p_2} dq
\int^{\min(q, 2p_1 - q)}_{\max(-q, q - 2p_2)} \frac{d\omega}{2\pi}
\exp(-p_2^2) \,
```
The statistical factor $\mathcal{F}[f]$ for this case is given by:
```{math}
:label:
\mathcal{F}[f] = f_2\;.
```

The integral over $\omega$ can be performed analytically using SageMath as follows:
:::{literalinclude} ../../notebooks/scalar_integrals.sage
:start-after: "# START"
:end-before: "# END"
:language: python
:caption: SageMath code to compute the integral over $\omega$.
:linenos:
:::
The first integral leads to:
```{math}
:label:
C[f_1]
= \frac{1}{2\nu_1}\frac{1}{{(2\pi)}^3}\frac{p_1}{16p_1^2}
= \frac{1}{2\nu_1}\frac{1}{{(2\pi)}^3}\frac{1}{16p_1}\;.
```

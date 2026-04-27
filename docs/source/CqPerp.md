(cqperp)=
# CqPerp

For computing $C(q_\perp)$, we start from the reduced collision integral {eq}`eq-reduced-collision-integral`.
```{math}
:label:
C[f_1] = & \frac{1}{2\nu_1} \int \frac{d^3p_2}{(2\pi)^3}
\frac{d^3 q}{(2\pi)^3}
\int_{-q}^{q} \frac{d\omega}{2\pi} \,
\frac{1}{16p_1 p_2 |\vec{p}_1 - \vec{q}| |\vec{p}_2 + \vec{q}|} \\
& \times (2\pi) \delta(p_1 - \omega - |\vec{p}_1 - \vec{q}|) \,
(2\pi) \delta(p_2 + \omega - |\vec{p}_2 + \vec{q}|) \,
|\mathcal{M}|^2 \, \mathcal{F}[f]\;.
```

Here we keep the matrix element $|\mathcal{M}|^2$ and statistical factor $\mathcal{F}[f]$ general. 
The statistical factor depends on the distribution functions:
```{math}
:label:
\mathcal{F}[f] = \mathcal{F}[f](\vec{p}_2, \vec{p}_2 + \vec{q})\;.
```

We now expand the first delta function for $\frac{q}{p_1} \to 0$
```{math}
:label:
\delta(p_1 - \omega - |\vec{p}_1 - \vec{q}|)
= \delta\left(\frac{\vec{p}_1 \cdot \vec{q}}{p_1} - \omega\right)
\;.
```
We define $q_\parallel = \frac{\vec{p}_1 \cdot \vec{q}}{p_1}$ such that the delta function becomes $\delta(q_\parallel - \omega)$.
The integration over $q$ is then expressed in terms of $q_\perp$ and $q_\parallel$ as
```{math}
:label:
\int d^3 q = \int d^2 q_\perp \, d q_\parallel
\;.
```
Using the delta function to perform the integral over $q_\parallel$ we get $q = \sqrt{q_\perp^2 + \omega^2}$ and $|\vec{p}_1 - \vec{q}| = p_1 - \omega$.
Now the bounds on $\omega$ are extended to infinity since we always have $-q \leq \omega \leq q$.

For the second delta function we use the same treatment as in {eq}`eq-jacobians`:
```{math}
:label:
\delta(p_2 + \omega - |\vec{p}_2 + \vec{q}|)
 = \frac{|\vec{p}_2 + \vec{q}|}{p_2 q} \,
  \delta\left(\cos\theta_{2q} - \left(\frac{\omega}{q} + \frac{\omega^2 - q^2}{2 p_2 q}\right)\right)
  \Theta\left(p_2 - \frac{q - \omega}{2}\right)\;,
```
Defining the $C(q_\perp)$
```{math}
:label:
C(q_\perp) = \frac{1}{{(2\pi)}^2} \frac{d\Gamma}{d^2 q_\perp}\;,
```
we obtain the general expression
```{math}
:label: eq-cqperp-general
C(q_\perp) = & \frac{1}{(2\pi)^2}
\int^{\infty}_{-\infty} \frac{d\omega}{2\pi} \,
\int^\infty_{\frac{q-\omega}{2}} dp_2 \, p_2^2
\int d\cos\theta_{2q} \, d\phi_{2q}
\\
& \times
\frac{|\mathcal{M}|^2}{16 p_1 (p_1 - \omega) p_2^2 q}
(2\pi) \delta\left(\cos\theta_{2q} - \left(\frac{\omega}{q} + \frac{\omega^2 - q^2}{2 p_2 q}\right)\right) \,
\mathcal{F}[f]\;.
```

Replacing $p_1 - \omega \approx p_1$ in the small-$q$ limit, the expression simplifies to:
```{math}
:label: eq-cqperp-simplified
C(q_\perp) = & \frac{1}{(2\pi)^2}
\int^{\infty}_{-\infty} \frac{d\omega}{2\pi} \,
\int^\infty_{\frac{q-\omega}{2}} dp_2 \,
\int d\cos\theta_{2q} \, d\phi_{2q}
\\
& \times
\frac{|\mathcal{M}|^2}{16 p_1^2 q}
(2\pi) \delta\left(\cos\theta_{2q} - \left(\frac{\omega}{q} + \frac{\omega^2 - q^2}{2 p_2 q}\right)\right) \,
\mathcal{F}[f]\;.
```

The matrix element $|\mathcal{M}|^2$ can depend on the kinematic variables $p_2$, $p_{2\parallel}$, $\omega$, $q_\perp$, and $q$.

## Basis Vectors

Similarly to {eq}`eq-basis-vectors-p1`, we introduce basis vectors for $\vec{q}$:
```{math}
:label:
\vec{e}^{\, 1}_{p_1} &= (\sin\theta_{p_1} \cos\phi_{p_1}, \sin\theta_{p_1} \sin\phi_{p_1}, \cos\theta_{p_1})\;, \\
\vec{e}^{\, 2}_{p_1} &= (\cos\theta_{p_1} \cos\phi_{p_1}, \cos\theta_{p_1} \sin\phi_{p_1}, -\sin\theta_{p_1})\;, \\
\vec{e}^{\,3}_{p_1} &= (-\sin\phi_{p_1}, \cos\phi_{p_1}, 0)\;.
```

Using $\cos\theta_{1q} = \frac{\omega}{q}$, we express the momentum $\vec{q}$ as:

```{math}
:label: eq-q-vector
\vec{q} = q \left[
 \cos\theta_{1q}\vec{e}^{\, 1}_{p_1}
 +\sin\theta_{1q} \left(\cos\phi_{1q}\vec{e}^{\, 2}_{p_1} + \sin\phi_{1q}\vec{e}^{\, 3}_{p_1}\right)
  \right]\;.
```

Similarly, we can express $\vec{p}_2$ in terms of the basis vectors defined by $\vec{q}$:

```{math}
:label: eq-basis-vectors-q
\vec{e}^{\, 1}_{q} &= (\sin\theta_{q} \cos\phi_{q}, \sin\theta_{q} \sin\phi_{q}, \cos\theta_{q})\;, \\
\vec{e}^{\, 2}_{q} &= (\cos\theta_{q} \cos\phi_{q}, \cos\theta_{q} \sin\phi_{q}, -\sin\theta_{q})\;, \\
\vec{e}^{\,3}_{q} &= (-\sin\phi_{q}, \cos\phi_{q}, 0)\;.
```

Leading to the expression:

```{math}
:label: eq-p2-vector
\vec{p}_2 = p_2 \left[
 \cos\theta_{2q}\vec{e}^{\, 1}_{q}
 +\sin\theta_{2q} \left(\cos\phi_{2q}\vec{e}^{\, 2}_{q} + \sin\phi_{2q}\vec{e}^{\, 3}_{q}\right)
  \right]\;.
```
The projection of $\vec{p}_2$ along $\vec{p}_1$ is simply
```{math}
:label:
p_{2\parallel} = (\vec{p}_2 \cdot \vec{e}^{\, 1}_{p_1})\;.
```

## Leading-Log Matrix Element

In the leading-log approximation for QCD, the matrix element takes the form:
```{math}
:label: eq-matrix-leading-log
|\mathcal{M}|^2_{\text{LL}} = 4 g^4 \frac{{\left(p_1 p_2 - \vec{p}_1 \cdot \vec{p}_2\right)}^2}{{\left(q^2 - \omega^2\right)}^2} 
= 4 g^4 \frac{{\left(p_2 - p_{2\parallel}\right)}^2}{q_\perp^4}\;,
```
where we used $p_1 p_2 - \vec{p}_1 \cdot \vec{p}_2 = p_1(p_2 - p_{2\parallel})$ and $q^2 - \omega^2 = q_\perp^2$.

For the statistical factor in QCD with gluons and quarks:
```{math}
:label: eq-stat-qcd
\mathcal{F}[f] = 
C_A f_g(\vec{p}_2) (1 + f_g(\vec{p}_2 + \vec{q}))
+ \frac{N_f}{2} f_q(\vec{p}_2) (1 - f_q(\vec{p}_2 + \vec{q}))
+ \frac{N_f}{2} f_{\bar{q}}(\vec{p}_2) (1 - f_{\bar{q}}(\vec{p}_2 + \vec{q}))
\;.
```

Substituting the leading-log matrix element into {eq}`eq-cqperp-simplified`:
```{math}
:label: eq-cqperp-leading-log
C(q_\perp)_{\text{LL}} = & \frac{1}{(2\pi)^2} \frac{1}{q_\perp^4}
\int^{\infty}_{-\infty} \frac{d\omega}{2\pi} \,
\int^\infty_{\frac{q-\omega}{2}} dp_2 \,
\int d\cos\theta_{2q} \, d\phi_{2q}
\\
& \times
\frac{g^4 {\left(p_2 - p_{2\parallel}\right)}^2}{4 q}
(2\pi) \delta\left(\cos\theta_{2q} - \left(\frac{\omega}{q} + \frac{\omega^2 - q^2}{2 p_2 q}\right)\right) \,
\mathcal{F}[f]\;.
```

Note that we have re-arranged the color factors into the statistical factor $\mathcal{F}[f]$ and stripped $C(q_\perp)$ of the Casimir factor $C_R$.

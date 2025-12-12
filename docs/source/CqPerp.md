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

Here the matrix element and statistical factors are simply 
```{math}
:label:
\frac{1}{2\nu_1}|\mathcal{M}|^2 &= 4 g^4 \frac{{\left(p_1 p_2 - \vec{p}_1 \cdot \vec{p}_2\right)}^2}{{\left(q^2 - \omega^2\right)}^2} \;,\\
\mathcal{F}[f] &= 
C_A f_g(\vec{p}_2) (1 + f_g(\vec{p}_2 + \vec{q}))
+ \frac{N_f}{2} f_q(\vec{p}_2) (1 - f_q(\vec{p}_2 + \vec{q}))
+ \frac{N_f}{2} f_{\bar{q}}(\vec{p}_2) (1 - f_{\bar{q}}(\vec{p}_2 + \vec{q}))
\;.
```
We have re-arranged the color factors into the statistical factor and stripped the $C(q_\perp)$ of the Casimir factor $C_R$.

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
 &= \frac{|\vec{p}_2 + \vec{q}|}{p_2 q} \,
  \delta\left(\cos\theta_{2q} - \left(\frac{\omega}{q} + \frac{\omega^2 - q^2}{2 p_2 q}\right)\right)
  \Theta\left(p_2 - \frac{q - \omega}{2}\right)\;,
```
Defining the $C(q_\perp)$
```{math}
:label:
C(q_\perp) = \frac{1}{{(2\pi)}^2} \frac{d\Gamma}{d^2 q_\perp}\;,
```
we obtain the expression
```{math}
:label:
C(q_\perp) = & \frac{1}{(2\pi)^2} \frac{1}{q_\perp^4}
\int^{\infty}_{-\infty} \frac{d\omega}{2\pi} \,
\int^\infty_{\frac{q-\omega}{2}} dp_2
d\cos\theta_{2q} d\phi_{2q}
\frac{4 g^4 {\left(p_1 p_2 - \vec{p}_1 \cdot \vec{p}_2\right)}^2p_2^2}{16p_1 (p_1 - \omega) p_2^2 q}
\\
& \times
(2\pi) \delta\left(\cos\theta_{2q} - \left(\frac{\omega}{q} + \frac{\omega^2 - q^2}{2 p_2 q}\right)\right) \,
\mathcal{F}[f]\;.
```
We replace $p_1 - \omega \approx p_1$ and simplify the expression to
```{math}
:label:
C(q_\perp) = & \frac{1}{(2\pi)^2} \frac{1}{q_\perp^4}
\int^{\infty}_{-\infty} \frac{d\omega}{2\pi} \,
\int^\infty_{\frac{q-\omega}{2}} dp_2
d\cos\theta_{2q} d\phi_{2q}
\frac{g^4 {\left(p_2 - p_{2\parallel}\right)}^2}{4q}
\\
& \times
(2\pi) \delta\left(\cos\theta_{2q} - \left(\frac{\omega}{q} + \frac{\omega^2 - q^2}{2 p_2 q}\right)\right) \,
\mathcal{F}[f]\;.
```

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


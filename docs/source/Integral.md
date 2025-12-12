(integral)=
# 2-to-2 Phase-Space Measure

The collision integral for $2 \leftrightarrow 2$ scattering processes is given by:

```{math}
:label: eq-collision-integral
C[f_1] = \frac{1}{2E_1} \int \frac{d^3p_2}{(2\pi)^3 2E_2} 
\frac{d^3p_3}{(2\pi)^3 2E_3} \frac{d^3p_4}{(2\pi)^3 2E_4} 
(2\pi)^4 \delta^{(4)}(p_1 + p_2 - p_3 - p_4) \, 
|\mathcal{M}|^2 \, \mathcal{F}[f]\;,
```

where the statistical factor is:

```{math}
:label: eq-statistical-factor
\mathcal{F}[f] = f_1 f_2 (1 \pm f_3)(1 \pm f_4) - f_3 f_4 (1 \pm f_1)(1 \pm f_2)
```

We first introduce the momentum transfer $q = p_1 - p_3$ and the energy transfer $\omega = E_1 - E_3$.
Writing the energy delta function as:

```{math}
:label: eq-energy-delta
\delta(E_1 + E_2 - E_3 - E_4) = \int d\omega \, \delta(E_1 - \omega - E_3) \, \delta(E_2 + \omega - E_4)
```

Using the delta functions we perform the integral over $p_4$ and shift the integral $p_3 \to q$.
The positivity of energies $E_3, E_4 \geq 0$ implies $-q \leq \omega \leq q$.
We obtain the reduced collision integral:

```{math}
:label: eq-reduced-collision-integral
C[f_1] = & \int \frac{d^3p_2}{(2\pi)^3}
\frac{d^3 q}{(2\pi)^3}
\int_{-q}^{q} \frac{d\omega}{2\pi} \,
\frac{1}{16p_1 p_2 |\vec{p}_1 - \vec{q}| |\vec{p}_2 + \vec{q}|} \\
& \times (2\pi) \delta(p_1 - \omega - |\vec{p}_1 - \vec{q}|) \,
(2\pi) \delta(p_2 + \omega - |\vec{p}_2 + \vec{q}|) \,
|\mathcal{M}|^2 \, \mathcal{F}[f]\;.
```

We use the angles $\cos\theta_{1q} = (\vec{p}_1 \cdot \vec{q})/(p_1 q)$ and $\cos\theta_{2q} = (\vec{p}_2 \cdot \vec{q})/(p_2 q)$ to perform the integrals over the delta functions.
The Jacobians are:

```{math}
:label: eq-jacobians
\delta(p_1 - \omega - |\vec{p}_1 - \vec{q}|)
 &= \frac{|\vec{p}_1 - \vec{q}|}{p_1 q} \,
  \delta\left(\cos\theta_{1q} - \left(\frac{\omega}{q} - \frac{\omega^2 - q^2}{2 p_1 q}\right)\right)
  \Theta\left(p_1 - \frac{q + \omega}{2}\right)\;,\\
\delta(p_2 + \omega - |\vec{p}_2 + \vec{q}|)
 &= \frac{|\vec{p}_2 + \vec{q}|}{p_2 q} \,
  \delta\left(\cos\theta_{2q} - \left(\frac{\omega}{q} + \frac{\omega^2 - q^2}{2 p_2 q}\right)\right)
  \Theta\left(p_2 - \frac{q - \omega}{2}\right)\;,
```

where the $\Theta$-functions together with the positivity of $p_1$ and $p_2$ ensures the conditions $|\cos\theta_{1q}| \leq 1$ and $|\cos\theta_{2q}| \leq 1$.

Combining the $\Theta$-functions, we obtain the integration limits:

```{math}
:label: eq-integration-limits
q &\leq p_1 + p_2\;, \\
\max(-q, q - 2p_2) &\leq \omega \leq \min(q, 2p_1 - q)\;.
```

Finally, we obtain the expression for the collision integral:

```{math}
:label: eq-final-collision-integral
C[f_1] = & \frac{1}{{(2\pi)}^6}
\int dp_2 \int_{0}^{p_1+p_2} dq
d\cos\theta_{2q} d\phi_{2q}
d\cos\theta_{1q} d\phi_{1q}
\int_{\min(q, 2p_1 - q)}^{\max(-q, q - 2p_2)} \frac{d\omega}{2\pi} \,
\frac{p_2^2 q^2}{16p_1^2 p_2^2 q^2} \\
& \times (2\pi) \delta\left(\cos\theta_{1q} - \left(\frac{\omega}{q} - \frac{\omega^2 - q^2}{2 p_1 q}\right)\right) \,
(2\pi) \delta\left(\cos\theta_{2q} - \left(\frac{\omega}{q} + \frac{\omega^2 - q^2}{2 p_2 q}\right)\right) \,
|\mathcal{M}|^2 \, \mathcal{F}[f]\;.
```

We introduce the basis vectors:

```{math}
:label: eq-basis-vectors-p1
\vec{e}^{\, 1}_{p_1} &= (\sin\theta_{p_1} \cos\phi_{p_1}, \sin\theta_{p_1} \sin\phi_{p_1}, \cos\theta_{p_1})\;, \\
\vec{e}^{\, 2}_{p_1} &= (\cos\theta_{p_1} \cos\phi_{p_1}, \cos\theta_{p_1} \sin\phi_{p_1}, -\sin\theta_{p_1})\;, \\
\vec{e}^{\,3}_{p_1} &= (-\sin\phi_{p_1}, \cos\phi_{p_1}, 0)\;.
```

Allowing us to express the momentum $\vec{q}$ as:

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

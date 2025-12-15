(process-traits)=
# Process Traits

The `ProcessTraits` system provides compile-time process selection for collision integrals. Each scattering process (e.g., $gg \to gg$, $qg \to qg$) is defined by a trait specialization that bundles together the relevant physics functions.

## Overview

A scattering process requires:

1. **Distribution functions** for particles 1-4
2. **Statistical factor** encoding quantum statistics (Bose-Einstein or Fermi-Dirac)
3. **Matrix element squared** $|\mathcal{M}|^2$

The traits system allows switching between processes at compile time with zero runtime overhead.

## Architecture

```
┌─────────────────┐     ┌──────────────────────────┐
│  Process Tag    │────▶│  ProcessTraits<Tag>      │
│  (gg_to_gg)     │     │  - dist1, dist2, ...     │
└─────────────────┘     │  - stat                  │
                        │  - matrix                │
                        └──────────────────────────┘
                                    │
                                    ▼
                        ┌──────────────────────────┐
                        │  CollisionIntegral<Tag>  │
                        └──────────────────────────┘
```

## Defining a Process

### Step 1: Create a Process Tag

A process tag is an empty struct that identifies the process:

```cpp
/// Tag for gg -> gg process
struct gg_to_gg {};

/// Tag for qg -> qg process  
struct qg_to_qg {};
```

### Step 2: Specialize ProcessTraits

Each process must specialize the `ProcessTraits` template using `constexpr auto`:

```cpp
template<>
struct ProcessTraits<gg_to_gg> {
    static constexpr auto dist1 = Distribution::gluon;
    static constexpr auto dist2 = Distribution::gluon;
    static constexpr auto dist3 = Distribution::gluon;
    static constexpr auto dist4 = Distribution::gluon;
    static constexpr auto stat = StatisticalFactor::BBBB;
    static constexpr auto matrix = MatrixElement::gg_gg;
};
```

The `constexpr auto` syntax lets the compiler deduce the function pointer type automatically, avoiding explicit type aliases.

### Step 3: Use in Collision Integral

```cpp
// Instantiate the integrand for gg -> gg
auto integrand = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>;

// Run the computation
IntegrateQCD::Compute<gg_to_gg>(integrand);
```

## Function Signatures

The traits reference functions with the following signatures:

### Distribution Function

```cpp
double distribution(double p, double cosTheta, double phi);
```

Computes the phase-space distribution $f(\vec{p})$ for a particle.

| Parameter | Description |
|-----------|-------------|
| `p` | Momentum magnitude $|\vec{p}|$ |
| `cosTheta` | Polar angle cosine $\cos\theta$ |
| `phi` | Azimuthal angle $\phi$ |

### Statistical Factor Function

```cpp
double statistical_factor(double f1, double f2, double f3, double f4);
```

Computes the statistical factor $\mathcal{F}[f]$ from {eq}`eq-statistical-factor`.

| Function | Formula | Particles |
|----------|---------|-----------|
| `BBBB` | $f_1 f_2 (1+f_3)(1+f_4) - f_3 f_4 (1+f_1)(1+f_2)$ | 4 bosons |
| `FBFB` | $f_1 f_2 (1-f_3)(1+f_4) - f_3 f_4 (1-f_1)(1+f_2)$ | fermion-boson |
| `FFFF` | $f_1 f_2 (1-f_3)(1-f_4) - f_3 f_4 (1-f_1)(1-f_2)$ | 4 fermions |

### Matrix Element Function

```cpp
double matrix_element(double s, double t, double u, 
                      double qt, double qu, 
                      double mDSqr, double mQSqr);
```

Computes the matrix element squared $|\mathcal{M}|^2$.

| Parameter | Description |
|-----------|-------------|
| `s, t, u` | Mandelstam variables |
| `qt` | Momentum transfer $|\vec{p}_1 - \vec{p}_3|$ |
| `qu` | Momentum transfer $|\vec{p}_2 - \vec{p}_3|$ |
| `mDSqr` | Debye mass squared $m_D^2$ |
| `mQSqr` | Quark thermal mass squared $m_q^2$ |

## Available Processes

### gg -> gg (Gluon-Gluon Scattering)

```cpp
struct gg_to_gg {};
```

All four particles are gluons with Bose-Einstein statistics.

```{math}
|\mathcal{M}|^2_{gg\to gg} = 4 g^4 d_A C_A^2 \left[
9 + \frac{(s-u)^2}{\bar{t}^2} + \frac{(s-t)^2}{\bar{u}^2} + \frac{(u-t)^2}{s^2}
\right]
```

where the screened propagators are:

```{math}
\bar{t} = t \frac{q_t^2 + \xi_0^2 m_D^2}{q_t^2}, \quad
\bar{u} = u \frac{q_u^2 + \xi_0^2 m_D^2}{q_u^2}
```

with screening parameter $\xi_0 = e^{5/6}/(2\sqrt{2})$.

### qg -> qg (Quark-Gluon Scattering)

```cpp
struct qg_to_qg {};
```

Particle 1 and 3 are quarks (fermions), particles 2 and 4 are gluons (bosons).

### qq -> qq (Quark-Quark Scattering)

```cpp
struct qq_to_qq {};
```

All four particles are quarks with Fermi-Dirac statistics.

## Adding a New Process

To add a new process (e.g., $q\bar{q} \to gg$):

1. **Define the tag:**
   ```cpp
   struct qqbar_to_gg {};
   ```

2. **Implement the matrix element** in `namespace MatrixElement`:
   ```cpp
   inline double qqbar_gg(double s, double t, double u,
                          double qt, double qu,
                          double mDSqr, double mQSqr) {
       // Implementation
   }
   ```

3. **Specialize the traits:**
   ```cpp
   template<>
   struct ProcessTraits<qqbar_to_gg> {
       static constexpr auto dist1 = Distribution::quark;
       static constexpr auto dist2 = Distribution::quark;  // antiquark
       static constexpr auto dist3 = Distribution::gluon;
       static constexpr auto dist4 = Distribution::gluon;
       static constexpr auto stat = StatisticalFactor::FFBB;
       static constexpr auto matrix = MatrixElement::qqbar_gg;
   };
   ```

4. **Use it:**
   ```cpp
   auto integrand = CollisionIntegralQCD::CollisionIntegral<qqbar_to_gg>;
   IntegrateQCD::Compute<qqbar_to_gg>(integrand);
   ```

## Design Rationale

### Why Compile-Time Selection?

1. **Zero overhead**: No virtual function calls or runtime dispatch
2. **Inlining**: Distribution and matrix element functions are inlined
3. **Type safety**: Invalid process configurations are caught at compile time
4. **Optimization**: Compiler can optimize for each specific process

### Why `constexpr auto`?

Using `constexpr auto` instead of explicit function pointer type aliases:

1. **Cleaner syntax**: No need to define `using DistFunc = double(*)(...)` etc.
2. **Type deduction**: Compiler automatically deduces the correct function pointer type
3. **Less boilerplate**: Adding new functions doesn't require updating type aliases
4. **C++17 feature**: Takes advantage of modern C++ capabilities

The `inline` functions in `Distribution::`, `StatisticalFactor::`, and `MatrixElement::` namespaces have their addresses as constant expressions, making them usable with `constexpr`.

## File Reference

| File | Contents |
|------|----------|
| `src/Process.hpp` | Process tags, traits, distributions, statistical factors, matrix elements |
| `src/Integral.cpp` | `CollisionIntegral<ProcessTag>` template |

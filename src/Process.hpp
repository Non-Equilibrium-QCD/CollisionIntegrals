#pragma once

#include <cmath>
#include "constants.cpp"

/**
 * @file Process.hpp
 * @brief Compile-time process definitions for collision integrals.
 * 
 * Each process defines:
 * - Distribution functions for particles 1-4
 * - Statistical factor
 * - Matrix element squared
 */

namespace Distribution {

/// Bose-Einstein distribution (equilibrium)
inline double bose(double p, double T) {
    return 1.0 / (std::exp(p / T) - 1.0);
}

/// Fermi-Dirac distribution (equilibrium)
inline double fermi(double p, double T) {
    return 1.0 / (std::exp(p / T) + 1.0);
}

/// Gluon distribution (non-equilibrium)
inline double gluon(double p, double cosTheta, double phi) {
    return std::exp(-p * p);// * std::exp(-cosTheta * cosTheta);
}

/// Quark distribution (non-equilibrium)
inline double quark(double p, double cosTheta, double phi) {
    return std::exp(-p * p);
}

} // namespace Distribution


namespace StatisticalFactor {

/// Boson-Boson -> Boson-Boson (e.g., gg -> gg)
inline double BBBB(double f1, double f2, double f3, double f4) {
    return f1 * f2 * (1.0 + f3) * (1.0 + f4) - f3 * f4 * (1.0 + f1) * (1.0 + f2);
}

/// Fermion-Boson -> Fermion-Boson (e.g., qg -> qg)
inline double FBFB(double f1, double f2, double f3, double f4) {
    return f1 * f2 * (1.0 - f3) * (1.0 + f4) - f3 * f4 * (1.0 - f1) * (1.0 + f2);
}

/// Fermion-Fermion -> Fermion-Fermion (e.g., qq -> qq)
inline double FFFF(double f1, double f2, double f3, double f4) {
    return f1 * f2 * (1.0 - f3) * (1.0 - f4) - f3 * f4 * (1.0 - f1) * (1.0 - f2);
}

} // namespace StatisticalFactor


namespace MatrixElement {

// SCREENING PARAMETER //
static const double xi0g = std::exp(5.0 / 6.0) / (2.0 * M_SQRT2);
static const double xi0q = std::exp(1.0) / M_SQRT2;


/// gg -> gg matrix element squared
inline double gg_gg(double s, double t, double u, 
                    double qt, double qu, 
                    double mDSqr, double mQSqr) {
    double tBar = t * (qt * qt + xi0g * xi0g * mDSqr) / (qt * qt);
    double uBar = u * (qu * qu + xi0g * xi0g * mDSqr) / (qu * qu);

    return 4.0 * g * g * g * g * dA * CA * CA * (
        9.0
        + (s - u) * (s - u) / (tBar * tBar)
        + (s - t) * (s - t) / (uBar * uBar)
        + (u - t) * (u - t) / (s * s)
    );
}

/// qg -> qg matrix element squared
inline double qg_qg(double s, double t, double u,
                    double qt, double qu,
                    double mDSqr, double mQSqr) {
    // TODO: implement
    return 0.0;
}

/// qq -> qq matrix element squared
inline double qq_qq(double s, double t, double u,
                    double qt, double qu,
                    double mDSqr, double mQSqr) {
    // TODO: implement
    return 0.0;
}

} // namespace MatrixElement


// ============================================================================
// Process Traits
// ============================================================================

/**
 * @brief Process trait template.
 * 
 * Specialize this for each scattering process to define:
 * - dist1, dist2, dist3, dist4: distribution functions
 * - stat: statistical factor
 * - matrix: matrix element squared
 */
template<typename ProcessTag>
struct ProcessTraits;

// ============================================================================
// Process Tags
// ============================================================================

/// Tag for gg -> gg process
struct gg_to_gg {};

/// Tag for qg -> qg process  
struct qg_to_qg {};

/// Tag for qq -> qq process
struct qq_to_qq {};

// ============================================================================
// Process Trait Specializations
// ============================================================================

template<>
struct ProcessTraits<gg_to_gg> {
    static constexpr auto dist1 = Distribution::gluon;
    static constexpr auto dist2 = Distribution::gluon;
    static constexpr auto dist3 = Distribution::gluon;
    static constexpr auto dist4 = Distribution::gluon;
    static constexpr auto stat = StatisticalFactor::BBBB;
    static constexpr auto matrix = MatrixElement::gg_gg;
    static constexpr double nuA = nuG;
};

template<>
struct ProcessTraits<qg_to_qg> {
    static constexpr auto dist1 = Distribution::quark;
    static constexpr auto dist2 = Distribution::gluon;
    static constexpr auto dist3 = Distribution::quark;
    static constexpr auto dist4 = Distribution::gluon;
    static constexpr auto stat = StatisticalFactor::FBFB;
    static constexpr auto matrix = MatrixElement::qg_qg;
    static constexpr double nuA = nuQ;
};

template<>
struct ProcessTraits<qq_to_qq> {
    static constexpr auto dist1 = Distribution::quark;
    static constexpr auto dist2 = Distribution::quark;
    static constexpr auto dist3 = Distribution::quark;
    static constexpr auto dist4 = Distribution::quark;
    static constexpr auto stat = StatisticalFactor::FFFF;
    static constexpr auto matrix = MatrixElement::qq_qq;
    static constexpr double nuA = nuQ;
};

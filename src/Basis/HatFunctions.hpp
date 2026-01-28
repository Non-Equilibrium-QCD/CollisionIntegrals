#pragma once

/**
 * @file HatFunctions.hpp
 * @brief Hat function basis for function expansions.
 *
 * HatFunctionBasis derives from BasisBase and implements piecewise linear
 * hat functions for basis expansions.
 */

#include <cmath>
#include <vector>

// ============================================================================
// Grid Parameters
// ============================================================================

struct GridParams {
    size_t Nx;
    double xmin, xmax;
};

/**
 * @brief Hat function basis implementation.
 *
 * Implements piecewise linear hat functions (tent functions) for each
 * grid dimension. Inherits grid data and utilities from BasisBase.
 */
struct HatFunctionBasis {
    // Per-derived static data (each Derived gets its own copy)
    static inline GridParams MomentumParams, CosThetaParams, PhiParams;
    static inline std::vector<double> MomentumVals, CosThetaVals, PhiVals;
    static inline std::vector<double> MomentumArea, CosThetaArea, PhiArea;

    /**
     * @brief Computes a hat function value at a given point.
     *
     * Hat functions are piecewise linear, with value 1 at the grid point
     * and 0 at neighboring grid points.
     *
     * @tparam Values Reference to the grid values vector
     * @param index The basis function index
     * @param variable The point at which to evaluate
     * @return The hat function value
     */
    template<std::vector<double>& Values>
    static inline double HatFct(size_t index, double variable) {
        size_t N = Values.size();
        size_t iL = (index > 0) ? index - 1 : 0;
        size_t iR = (index < N - 1) ? index + 1 : N - 1;

        if (variable <= Values[iL] || variable >= Values[iR]) {
            return 0.0;
        }

        if (variable < Values[index]) {
            return (variable - Values[iL]) / (Values[index] - Values[iL]);
        }

        return (Values[iR] - variable) / (Values[iR] - Values[index]);
    }

    /**
     * @brief Computes a hat function Area at a given point.
     *
     * \f[
     *     A[i] = \int_{x_{i-1}}^{x_{i+1}}~ hat_i(x) ~dx = \frac{(x_{i+1} - x_{i-1})}{2}
     * \f]
     *
     * @tparam Values Reference to the grid values vector
     * @param index The basis function index
     * @return The hat function area
     */
    template<std::vector<double>& Values>
    static inline double HatArea(size_t index) {
        size_t N = Values.size();
        if (N == 1) {
            return 1.0;
        }

        if (index == 0) {
            return (Values[index + 1] - Values[index]) / 2.0;
        }

        if (index == N - 1) {
            return (Values[index] - Values[index - 1]) / 2.0;
        }

        return (Values[index + 1] - Values[index - 1]) / 2.0;
    }

    /**
     * @brief Computes a hat function Area for p*f constant at a given point.
     *
     * \f[
     *     A[i] = \int_{x_{i-1}}^{x_{i+1}}~ x * hat_i(x) ~dx
     * \f]
     *
     * @tparam Values Reference to the grid values vector
     * @param index The basis function index
     * @return The hat function area
     */
    template<std::vector<double>& Values>
    static inline double HatArea_x(size_t index) {
        size_t N = Values.size();

        if (index == 0) {
            double xC = Values[index];
            double xR = Values[index + 1];
            return 1.0 / 6.0 * (-2.0 * xC * xC + xR * xC + xR * xR);
        }

        if (index == N - 1) {
            double xL = Values[index - 1];
            double xC = Values[index];
            return 1.0 / 6.0 * (2.0 * xC * xC - xL * xC - xL * xL);
        }

        double xL = Values[index - 1];
        double xC = Values[index];
        double xR = Values[index + 1];

        return 1.0 / 6.0 * (xL + xC + xR) * (xR - xL);
    }

    static inline double MomentumBasis(size_t index, double p) {
        return HatFct<MomentumVals>(index, p);
    }

    static inline double CosThetaBasis(size_t index, double cosTheta) {
        return HatFct<CosThetaVals>(index, cosTheta);
    }

    static inline double PhiBasis(size_t index, double phi) {
        return HatFct<PhiVals>(index, phi);
    }

    /**
     * @brief Generates a linear grid of values.
     * @param params Grid parameters (Nx, xmin, xmax)
     * @param values Output vector of grid values
     */
    static inline void LinearGrid(GridParams params,
                                  std::vector<double> &values) {
        values.resize(params.Nx);
        double step = (params.xmax - params.xmin) / (params.Nx - 1);

        for (size_t i = 0; i < params.Nx; i++) {
            values[i] = params.xmin + i * step;
        }
    }

    /**
     * @brief Finds the grid index for a given value (inverse of LinearGrid).
     * @param params Grid parameters
     * @param value The value to find index for
     * @return The grid index
     */
    static inline size_t InverseLinearGrid(GridParams params, double value) {
        double step = (params.xmax - params.xmin) / (params.Nx - 1);
        return static_cast<size_t>(std::floor((value - params.xmin) / step));
    }

    /**
     * @brief Sets up the basis grids.
     * @param mp Momentum grid parameters
     * @param cp CosTheta grid parameters
     * @param pp Phi grid parameters
     */
    static inline void Setup(GridParams mp, GridParams cp, GridParams pp) {
        MomentumParams = mp;
        CosThetaParams = cp;
        PhiParams = pp;
        LinearGrid(MomentumParams, MomentumVals);
        LinearGrid(CosThetaParams, CosThetaVals);
        LinearGrid(PhiParams, PhiVals);


        // Compute Areas
        MomentumArea.resize(MomentumParams.Nx);
        CosThetaArea.resize(CosThetaParams.Nx);
        PhiArea.resize(PhiParams.Nx);

        for (size_t i = 0; i < MomentumParams.Nx; i++) {
            double p = MomentumVals[i];
            double Measure = p / std::pow(2.0 * M_PI, 3);
            MomentumArea[i] = HatArea_x<MomentumVals>(i) * Measure;
        }

        for (size_t i = 0; i < CosThetaParams.Nx; i++) {
            CosThetaArea[i] = HatArea<CosThetaVals>(i);
        }

        for (size_t i = 0; i < PhiParams.Nx; i++) {
            PhiArea[i] = HatArea<PhiVals>(i);
        }
    }
};

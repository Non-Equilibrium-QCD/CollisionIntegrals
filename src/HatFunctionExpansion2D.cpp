#pragma once
#include <array>
#include <cmath>
#include "Basis/HatFunctions.hpp"
#include "ProcessTraits.hpp"
#include "constants.cpp"
#include <gsl/gsl_monte_vegas.h>
#include <vector>
#include "GSLVEGAS.cpp"
#include "GAUSSQUADRATURE.cpp"
#include <omp.h>
#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/color.h>

typedef QUAD<3, 52> GaussQuad3D;

namespace {
struct GSLARGS {
    size_t indexP, indexCosTheta;
    double pMin, pMax;
    double cosThetaMin, cosThetaMax;
    double (*fct)(double *, const GSLARGS &);
};
}

namespace FunctionExpansion {
constexpr size_t dimensions = 3;


/**
 * @brief Computes the expansion of a function f(p, cosTheta, phi)
 *        in terms of basis functions.
 * @tparam Basis Basis class derived from BasisBase, must define:
 *         - MomentumBasis(index, p) : basis functions in momentum space
 *         - CosThetaBasis(index, cosTheta) : basis functions in angular space
 *         - PhiBasis(index, phi) : basis functions in azimuthal space
 * @tparam fct The function to expand
 * @param x[3] Array of variables (p, cosTheta, phi)
 * @param args Struct containing basis indices
 * @return The integrand value for the expansion coefficient
 */
template<double (*fct)(double, double, double)>
inline double Expansion(std::array<double, 3> x, const GSLARGS& args) {
    size_t indexP = args.indexP;
    size_t indexCosTheta = args.indexCosTheta;
    double pMin = args.pMin;
    double pMax = args.pMax;
    double cosThetaMin = args.cosThetaMin;
    double cosThetaMax = args.cosThetaMax;
    double phiMin = 0.0;
    double phiMax = 2.0 * M_PI;
    double p = pMin + (pMax - pMin) * x[0]; // Map [0,1] to [pMin, pMax]
    double cosTheta = cosThetaMin + (cosThetaMax - cosThetaMin) *
                      x[1]; // Map [0,1] to [cosThetaMin, cosThetaMax]
    double phi = phiMin + (phiMax - phiMin) * x[2]; // Map [0,1] to [phiMin, phiMax]
    double Jacobian = (pMax - pMin) * (cosThetaMax - cosThetaMin) *
                      (phiMax - phiMin);
    double Measure = p * p / std::pow(2.0 * M_PI, 3);

    return HatFunctionBasis::MomentumBasis(indexP, p) *
           HatFunctionBasis::CosThetaBasis(indexCosTheta, cosTheta) *
           fct(p, cosTheta, phi) * Jacobian * Measure;
}

/**
 * @brief Computes the expansion coefficients of a function.
 * @tparam Basis class derived from BasisBase
 * @tparam Fct The function to expand: double(double p, double cosTheta, double phi)
 * @param OutputFile Path to output file
 */
template<double (*Fct)(double, double, double)>
std::vector<double> Compute() {
    auto Integrand = FunctionExpansion::Expansion<Fct>;

    size_t Np = HatFunctionBasis::MomentumParams.Nx;
    size_t Ncos = HatFunctionBasis::CosThetaParams.Nx;
    std::vector<double> Results(Np * Ncos);

    for (size_t i = 0; i < Np; i++) {
        double pMin = HatFunctionBasis::MomentumVals[
                          (i == 0) ? 0 : i - 1];
        double pMax = HatFunctionBasis::MomentumVals[
                          (i == Np - 1) ? Np - 1 : i + 1];

        for (size_t j = 0; j < Ncos; j++) {
            double cosThetaMin = HatFunctionBasis::CosThetaVals[
                                     (j == 0) ? 0 : j - 1];
            double cosThetaMax = HatFunctionBasis::CosThetaVals[
                                     (j == Ncos - 1) ? Ncos - 1 : j + 1];

            GSLARGS args;
            args.indexP = i;
            args.indexCosTheta = j;
            args.pMin = pMin;
            args.pMax = pMax;
            args.cosThetaMin = cosThetaMin;
            args.cosThetaMax = cosThetaMax;
            double result = GaussQuad3D::Integrate(Integrand, args);
            Results[j + i * Ncos] = result;

            if (i % (Np / 4) == 0 && j == 0) {
                fmt::println(stderr,
                             "Progress: {}/{} p-values computed \r",
                             i, Np);
            }
        }
    }

    return Results;
}

}

namespace CollisionExpansion {
std::vector<GSLVEGAS> vegasIntegrators;
constexpr size_t dimensions = 8;

template<typename ProcessTag>
inline double CollisionIntegral(double *x, const GSLARGS& args) {
    using Traits = ProcessTraits<ProcessTag>;

    size_t indexP = args.indexP;
    size_t indexCosTheta = args.indexCosTheta;
    // double pMin = args.pMin;
    // double pMax = args.pMax;
    // double cosThetaMin = args.cosThetaMin;
    // double cosThetaMax = args.cosThetaMax;
    // double phiMin = args.phiMin;
    // double phiMax = args.phiMax;
    // double p = pMin + (pMax - pMin) * x[0]; // Map [0,1] to [pMin, pMax]
    // double cosTheta = cosThetaMin + (cosThetaMax - cosThetaMin) * x[1]; // Map [0,1] to [cosThetaMin, cosThetaMax]
    // double phi = phiMin + (phiMax - phiMin) * x[2]; // Map [0,1] to [phiMin, phiMax]
    // double Jacobian = (pMax - pMin) * (cosThetaMax - cosThetaMin) * (phiMax - phiMin);

    // GET p1 ANGULAR COORDINATES //
    // SAMPLE p1 //
    double pMax = args.pMax;
    double pMin = args.pMin;
    double cosThetaMax = args.cosThetaMax;
    double cosThetaMin = args.cosThetaMin;
    double p1 = pMin + (pMax - pMin) * x[0]; // Map [0,1] to [pMin, pMax]
    double cosTheta1 = cosThetaMin + (cosThetaMax - cosThetaMin) * x[1]; // Map [0,1] to [cosThetaMin, cosThetaMax]
    double Sin1 = sqrt(1.0 - cosTheta1 * cosTheta1);
    double phi1 = 2.0 * M_PI * x[2];
    double Jacobian = (2.0 * M_PI) * (pMax - pMin) * (cosThetaMax - cosThetaMin);

    double p1x = p1 * Sin1 * cos(phi1);
    double p1y = p1 * Sin1 * sin(phi1);
    double p1z = p1 * cosTheta1;

    // SET BASIS VECTORS WRT p1 //
    double ep1[3]; double ep2[3]; double ep3[3];

    ep1[0] = Sin1 * cos(phi1);
    ep1[1] = Sin1 * sin(phi1);
    ep1[2] = cosTheta1;

    ep2[0] = cosTheta1 * cos(phi1);
    ep2[1] = cosTheta1 * sin(phi1);
    ep2[2] = -Sin1;

    ep3[0] = -sin(phi1);
    ep3[1] = cos(phi1);
    ep3[2] = 0.0;

    // SAMPLE p2 //
    double p2 = x[3] / (1.0 - x[3]) ; // maps [0,1] to [0,inf)

    if (x[3] == 1.0) {
        return 0.0;
    }

    // x[0] = p2 / (1.0 + p2);
    Jacobian *= (1.0 + p2) * (1.0 + p2);

    // SAMPLE q //
    double qMin = 0.0; double qMax = (p1 + p2);
    double q = qMin + (qMax - qMin) * x[4];
    Jacobian *= (qMax - qMin);

    // SAMPLE w //
    double wMin = std::max(-q, q - 2.0 * p2);
    double wMax = std::min(+q, 2.0 * p1 - q);
    double w = wMin + (wMax - wMin) * x[5];
    Jacobian *= (wMax - wMin);

    // SAMPLE Phi1q //
    double Phi1Q = 2.0 * M_PI * x[6];
    Jacobian *= (2.0 * M_PI);

    // SAMPLE Phi2q //
    double Phi2Q = 2.0 * M_PI * x[7];
    Jacobian *= (2.0 * M_PI);

    // GET Cos1Q AND Cos2Q FROM ENERGY CONSERVATION CONSTRAINTS //
    // double Cos1Q = -((p1 - w) * (p1 - w) - p1 * p1 - q * q) / (2.0 * p1 * q);
    double Cos1Q = w / q - (w * w - q * q) / (2.0 * p1 * q);
    double Sin1Q = sqrt(1.0 - Cos1Q * Cos1Q);
    // double Cos2Q = +((p2 + w) * (p2 + w) - p2 * p2 - q * q) / (2.0 * p2 * q);
    double Cos2Q = w / q + (w * w - q * q) / (2.0 * p2 * q);
    double Sin2Q = sqrt(1.0 - Cos2Q * Cos2Q);

    // SET q VECTOR //
    double qx = q * (Cos1Q * ep1[0] + Sin1Q * cos(Phi1Q) * ep2[0] + Sin1Q * sin(
                         Phi1Q) * ep3[0]);
    double qy = q * (Cos1Q * ep1[1] + Sin1Q * cos(Phi1Q) * ep2[1] + Sin1Q * sin(
                         Phi1Q) * ep3[1]);
    double qz = q * (Cos1Q * ep1[2] + Sin1Q * cos(Phi1Q) * ep2[2] + Sin1Q * sin(
                         Phi1Q) * ep3[2]);

    // DETERMINE ANGLES //
    double CosQ = qz / q; double SinQ = sqrt(1.0 - CosQ * CosQ);
    double PhiQ = atan2(qy, qx);

    // SET BASIS VECTORS WRT q //
    double eq1[3]; double eq2[3]; double eq3[3];

    eq1[0] = SinQ * cos(PhiQ);
    eq1[1] = SinQ * sin(PhiQ);
    eq1[2] = CosQ;

    eq2[0] = CosQ * cos(PhiQ);
    eq2[1] = CosQ * sin(PhiQ);
    eq2[2] = -SinQ;

    eq3[0] = -sin(PhiQ);
    eq3[1] = cos(PhiQ);
    eq3[2] = 0.0;

    // SET p2 VECTOR //
    double p2x = p2 * (Cos2Q * eq1[0] + Sin2Q * cos(Phi2Q) * eq2[0] + Sin2Q * sin(
                           Phi2Q) * eq3[0]);
    double p2y = p2 * (Cos2Q * eq1[1] + Sin2Q * cos(Phi2Q) * eq2[1] + Sin2Q * sin(
                           Phi2Q) * eq3[1]);
    double p2z = p2 * (Cos2Q * eq1[2] + Sin2Q * cos(Phi2Q) * eq2[2] + Sin2Q * sin(
                           Phi2Q) * eq3[2]);

    // SET p3 VECTOR //
    double p3x = p1x - qx;
    double p3y = p1y - qy;
    double p3z = p1z - qz;

    double p3 = sqrt(p3x * p3x + p3y * p3y + p3z * p3z);

    // SET p4 VECTOR //
    double p4x = p2x + qx;
    double p4y = p2y + qy;
    double p4z = p2z + qz;

    double p4 = sqrt(p4x * p4x + p4y * p4y + p4z * p4z);
    double xMax = HatFunctionBasis::MomentumParams.xmax;

    if (p2 >= xMax || p3 >= xMax || p4 >= xMax) {
        return 0.0;
    }

    // GET MANDELSTAM VARIABLES //
    double s = +2.0 * (p1 * p2 - (p1x * p2x + p1y * p2y + p1z * p2z));
    double t = w * w - q * q;
    double u = -s - t;

    // SET qt AND qu //
    double q13 = q;
    double q23 = std::sqrt((p2x - p3x) * (p2x - p3x) + (p2y - p3y) * (p2y - p3y) +
                           (p2z - p3z) * (p2z - p3z));

    double cosTheta2 = p2z / p2;
    double phi2 = atan2(p2y, p2x);
    double cosTheta3 = p3z / p3;
    double phi3 = atan2(p3y, p3x);
    double cosTheta4 = p4z / p4;
    double phi4 = atan2(p4y, p4x);
    double f1 = Traits::dist1(p1, cosTheta1, phi1);
    double f2 = Traits::dist2(p2, cosTheta2, phi2);
    double f3 = Traits::dist3(p3, cosTheta3, phi3);
    double f4 = Traits::dist4(p4, cosTheta4, phi4);

    // 1 + 2 -> 3 + 4
    double s1 = Traits::stat(f1, f2, f3, f4);
    double M1 = Traits::matrix(s, t, u, q13, q23, mDSqr, mQSqr);

    // // 2 + 1 -> 3 + 4
    // double s2 = Traits::stat(f2, f1, f3, f4);
    // double M2 = Traits::matrix(s, u, t, q23, q13, mDSqr, mQSqr);
    //
    // // 3 + 4 -> 1 + 2
    // double s3 = Traits::stat(f3, f4, f1, f2);
    // double M3 = Traits::matrix(s, t, u, q13, q23, mDSqr, mQSqr);
    // // 4 + 3 -> 1 + 2
    // double s4 = Traits::stat(f4, f3, f1, f2);
    // double M4 = Traits::matrix(s, u, t, q23, q13, mDSqr, mQSqr);

    double j  = Jacobian / (16.0 * p1 * p1);
    constexpr double TwoPi = power_recursive<double, 5>(2.0 * M_PI);
    constexpr double c  = 1.0 / (2.0 * Traits::nuA * TwoPi);

    if (std::abs(u) < 1e-12 || std::abs(t) < 1e-12) {
        return 0.0;
    }

    if (!std::isfinite(s1) || !std::isfinite(M1) || !std::isfinite(j)
            || !std::isfinite(c)) {
        fmt::println(stderr, "NaN encountered in integrand:"
                             " s1={}, M1={}, j={}, c={}",
                     s1, M1, j, c);
        return 0.0;
    }


    double sM = 0.0;
    sM += HatFunctionBasis::MomentumBasis(indexP, p1) *
          HatFunctionBasis::CosThetaBasis(indexCosTheta, cosTheta1) *
          s1 * M1;
    // sM += HatFunctionBasis::MomentumBasis(indexP, p2) *
    //       HatFunctionBasis::CosThetaBasis(indexCosTheta, cosTheta2) *
    //       s2 * M2;
    // sM -= HatFunctionBasis::MomentumBasis(indexP, p3) *
    //       HatFunctionBasis::CosThetaBasis(indexCosTheta, cosTheta3) *
    //       s3 * M3;
    // sM -= HatFunctionBasis::MomentumBasis(indexP, p4) *
    //       HatFunctionBasis::CosThetaBasis(indexCosTheta, cosTheta4) *
    //       s4 * M4;

    double Measure = p1 * p1 / std::pow(2.0 * M_PI, 3);
    return j * c * sM * Measure;
}

template<typename ProcessTag>
void Compute(std::string OutputFile, int Ncalls) {
    using Traits = ProcessTraits<ProcessTag>;
    auto Integrand = CollisionIntegral<ProcessTag>;

    size_t Np = HatFunctionBasis::MomentumParams.Nx;
    size_t Ncos = HatFunctionBasis::CosThetaParams.Nx;
    std::vector<double> Results(Np * Ncos), Errors(Np * Ncos);

    #pragma omp parallel for collapse(2)

    for (size_t i = 0; i < Np; i++) {
        for (size_t j = 0; j < Ncos; j++) {
            GSLARGS args;
            args.indexP = i;
            args.indexCosTheta = j;
            args.pMin = HatFunctionBasis::MomentumVals[
                            (i == 0) ? 0 : i - 1];;
            args.pMax = HatFunctionBasis::MomentumVals[
                            (i == Np - 1) ? Np - 1 : i + 1];
            args.cosThetaMin = HatFunctionBasis::CosThetaVals[
                                   (j == 0) ? 0 : j - 1];;
            args.cosThetaMax = HatFunctionBasis::CosThetaVals[
                                   (j == Ncos - 1) ? Ncos - 1 : j + 1];
            args.fct = Integrand;
            size_t tID = omp_get_thread_num();
            double result, error;
            vegasIntegrators[tID].integrate(args, Ncalls, &result, &error);
            Results[j + i * Ncos] = result;
            Errors[j + i * Ncos] = error;
            #pragma omp critical
            {
                if (i % (Np / 4) == 0 && j == 0) {
                    fmt::println(stderr,
                                 "Progress: {}/{} p-values computed \r",
                                 i, Np);
                }
            }
        }
    }

    // OUTPUT RESULTS //
    auto file = fmt::output_file(OutputFile);
    file.print(
        "# Each double block corresponds to a fixed Phi value.\n"
        "# Each block corresponds to a fixed p value. \n"
        "# Within each block, rows correspond to cosTheta values.\n"
        "# Gnuplot can plot this using index: splot 'filename' u 1:2:3 index 0 \n"
        "# for first Phi block, index 1 for second, etc.\n"
        "# in python you can use data = np.loadtxt('filename') to read the data.\n"
        "# Then reshape using data.reshape((Nphi, Np, Ncos, 3)) \n"
        "# the last 3 corresponds to (p, cosTheta, ExpansionCoefficient)\n"
    );

    auto dist1 = FunctionExpansion::Compute<Traits::dist1>();
    std::vector<double> dist2, dist3, dist4;

    // COMPTIME CHECK IF DISTRIBUTIONS ARE THE SAME TO AVOID REDUNDANT CALCULATIONS //
    if constexpr(Traits::dist1 != Traits::dist2) {
        dist2 = FunctionExpansion::Compute<Traits::dist2>();

    } else {
        dist2 = dist1;
    }

    if constexpr(Traits::dist1 != Traits::dist3 && Traits::dist2 != Traits::dist3) {
        dist3 = FunctionExpansion::Compute<Traits::dist3>();

    } else if constexpr(Traits::dist1 == Traits::dist3) {
        dist3 = dist1;

    } else {
        dist3 = dist2;
    }

    if constexpr(Traits::dist1 != Traits::dist4 && Traits::dist2 != Traits::dist4
                 && Traits::dist3 != Traits::dist4) {
        dist4 = FunctionExpansion::Compute<Traits::dist4>();

    } else if constexpr(Traits::dist1 == Traits::dist4) {
        dist4 = dist1;

    } else if constexpr(Traits::dist2 == Traits::dist4) {
        dist4 = dist2;

    } else {
        dist4 = dist3;
    }

    file.print("# p cosTheta n1_i n2_i n3_i n4_i C_i C_iErr Area\n");

    for (size_t i = 0; i < Np; i++) {
        for (size_t j = 0; j < Ncos; j++) {
            file.print("{} {} {} {} {} {} {} {} {}\n",
                       HatFunctionBasis::MomentumVals[i],
                       HatFunctionBasis::CosThetaVals[j],
                       dist1[j + i * Ncos],
                       dist2[j + i * Ncos],
                       dist3[j + i * Ncos],
                       dist4[j + i * Ncos],
                       Results[j + i * Ncos],
                       Errors[j + i * Ncos],
                       (2.0 * M_PI) *
                       HatFunctionBasis::MomentumArea[i] *
                       HatFunctionBasis::CosThetaArea[j]);
        }

        file.print("\n");
    }


}

void Setup() {
    vegasIntegrators.reserve(omp_get_max_threads());

    for (int i = 0; i < omp_get_max_threads(); i++) {
        vegasIntegrators.push_back(GSLVEGAS(CollisionExpansion::dimensions));
    }
}
}

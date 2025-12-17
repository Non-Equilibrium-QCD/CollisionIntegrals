#include "../src/Integral.cpp"

namespace DistributionQCD {
/// Gluon distribution (non-equilibrium)
inline double gluon(double p, double cosTheta, double phi) {
    (void)cosTheta; (void)phi;
    return std::exp(-p * p);// * std::exp(-cosTheta * cosTheta);
}
}

namespace StatisticalQCD {
inline double BBBB(double f1, double f2, double f3, double f4) {
    return f1 * f2 * (1.0 + f3) * (1.0 + f4) - f3 * f4 * (1.0 + f1) * (1.0 + f2);
}
}


namespace MatrixElementQCD {

// SCREENING PARAMETER //
static const double xi0g = std::exp(5.0 / 6.0) / (2.0 * M_SQRT2);
static const double xi0q = std::exp(1.0) / M_SQRT2;


/// gg -> gg matrix element squared
inline double gg_gg(double s, double t, double u,
                    double qt, double qu,
                    double mDSqr, double mQSqr) {
    (void)mQSqr; // not used in gg->gg
    double tBar = t * (qt * qt + xi0g * xi0g * mDSqr) / (qt * qt);
    double uBar = u * (qu * qu + xi0g * xi0g * mDSqr) / (qu * qu);

    return 4.0 * g * g * g * g * dA * CA * CA * (
               9.0
               + (s - u) * (s - u) / (tBar * tBar)
               + (s - t) * (s - t) / (uBar * uBar)
               + (u - t) * (u - t) / (s * s)
           );
}
}

/// Tag for gg -> gg process
struct QCDgg_gg {};

template<>
struct ProcessTraits<QCDgg_gg> {
    static constexpr auto dist1 = DistributionQCD::gluon;
    static constexpr auto dist2 = DistributionQCD::gluon;
    static constexpr auto dist3 = DistributionQCD::gluon;
    static constexpr auto dist4 = DistributionQCD::gluon;
    static constexpr auto stat  = StatisticalQCD::BBBB;
    static constexpr auto matrix = MatrixElement::gg_gg;
    static constexpr double nuA = nuG;
};

int main() {
    IntegrateQCD::Setup();
    IntegrateQCD::Compute<QCDgg_gg>("OUTPUT/QCDgg_gg.dat");
    return EXIT_SUCCESS;
}

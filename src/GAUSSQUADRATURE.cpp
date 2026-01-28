#pragma once

#include <array>
#include <omp.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>
template<unsigned Dim, unsigned NperDim>
struct QUAD_STATIC {

    // For odd NperDim, 0 is also a point
    static constexpr unsigned Half = (NperDim + 1) / 2;

    // Number of unique symmetric points
    using boost_type =
        boost::math::quadrature::detail::gauss_detail<double, NperDim, Half>;

    static inline std::array<double, NperDim> make_points(const
            std::vector<double> &bst, bool negative = false) {
        std::array<double, NperDim> pts{};
        unsigned non_zero_start = 1;

        if (NperDim & 1u) {
            pts[0] = bst[0];

        } else {
            non_zero_start = 0;
        }

        for (unsigned i = non_zero_start; i < Half; ++i) {
            pts[i] = bst[i];
            unsigned mirror_index = (NperDim & 1u) ? (NperDim - i) : (NperDim - 1 - i);
            pts[mirror_index] = negative ? -bst[i] : bst[i];
        }

        return pts;
    }

    static inline std::array<double, NperDim> points = make_points(
            boost_type::abscissa(), true);
    static inline std::array<double, NperDim> weights = make_points(
            boost_type::weights());
};


template<unsigned Dim, unsigned NperDim>
struct QUAD : QUAD_STATIC<Dim, NperDim> {
    using QUAD_STATIC<Dim, NperDim>::points;
    using QUAD_STATIC<Dim, NperDim>::weights;

    // Recursive loop nesting template using if constexpr
    // CurDim: current dimension being iterated (0 to Dim-1)
    template<unsigned CurDim, typename F>
    inline static void Integrate_recursive(
        const F& f,
        std::array<double, Dim> &x,
        double &sum,
        double wprod)
    {
        if constexpr(CurDim == Dim) {
            // Base case: all coordinates are set, evaluate function
            sum += wprod * f(x);

        } else {
            // Recursive case: loop over current dimension
            for (unsigned i = 0; i < NperDim; ++i) {
                double x_val = 0.5 + 0.5 * points[i];
                double w_val = weights[i] * 0.5;
                x[CurDim] = x_val;

                Integrate_recursive < CurDim + 1, F > (
                    f, x, sum, wprod * w_val
                );
            }
        }
    }

    template<typename F, typename... Args>
    static double Integrate(const F& f, Args... args) {
        double sum = 0.0;
        std::array<double, Dim> x;

        #pragma omp parallel for reduction(+:sum)

        for (unsigned i0 = 0; i0 < NperDim; ++i0) {
            double x_val = 0.5 + 0.5 * points[i0];
            double w_val = weights[i0] * 0.5;
            std::array<double, Dim> x_local = x;
            x_local[0] = x_val;

            Integrate_recursive<1, F>(f, x_local, sum, w_val);
        }

        return sum;
    }
};

// SPECIALIZATION for 1D
template<unsigned NperDim>
struct QUAD<1, NperDim> : QUAD_STATIC<1, NperDim> {
    using QUAD_STATIC<1, NperDim>::points;
    using QUAD_STATIC<1, NperDim>::weights;
    template<typename F, typename... Args>
    static double Integrate(const F& f, Args... args) {
        double sum = 0.0;

        for (unsigned i0 = 0; i0 < NperDim; ++i0) {
            std::array<double, 1> x = {0.5 + 0.5 * points[i0]};
            sum += weights[i0] * 0.5 * f(x, args...);
        }

        return sum;
    }
};

// SPECIALIZATION for 2D
template<unsigned NperDim>
struct QUAD<2, NperDim> : QUAD_STATIC<2, NperDim> {
    using QUAD_STATIC<2, NperDim>::points;
    using QUAD_STATIC<2, NperDim>::weights;
    template<typename F, typename... Args>
    static double Integrate(const F& f, Args... args) {
        double sum = 0.0;
        #pragma omp parallel for reduction(+:sum)

        for (unsigned i0 = 0; i0 < NperDim; ++i0) {
            double w0 = weights[i0] * 0.5;
            double x0 = 0.5 + 0.5 * points[i0];

            for (unsigned i1 = 0; i1 < NperDim; ++i1) {
                std::array<double, 2> x = {x0, 0.5 + 0.5 * points[i1]};
                sum += w0 * weights[i1] * 0.5 * f(x, args...);
            }
        }

        return sum;
    }
};

// SPECIALIZATION for 3D
template<unsigned NperDim>
struct QUAD<3, NperDim> : QUAD_STATIC<3, NperDim> {
    using QUAD_STATIC<3, NperDim>::points;
    using QUAD_STATIC<3, NperDim>::weights;
    template<typename F, typename... Args>
    static double Integrate(const F& f, Args... args) {
        double sum = 0.0;
        #pragma omp parallel for reduction(+:sum)

        for (unsigned i0 = 0; i0 < NperDim; ++i0) {
            double w0 = weights[i0] * 0.5;
            double x0 = 0.5 + 0.5 * points[i0];

            for (unsigned i1 = 0; i1 < NperDim; ++i1) {
                double w01 = w0 * weights[i1] * 0.5;
                double x1 = 0.5 + 0.5 * points[i1];

                for (unsigned i2 = 0; i2 < NperDim; ++i2) {
                    std::array<double, 3> x = {x0, x1, 0.5 + 0.5 * points[i2]};
                    sum += w01 * weights[i2] * 0.5 * f(x, args...);
                }
            }
        }

        return sum;
    }
};

// SPECIALIZATION for 4D
template<unsigned NperDim>
struct QUAD<4, NperDim> : QUAD_STATIC<4, NperDim> {
    using QUAD_STATIC<4, NperDim>::points;
    using QUAD_STATIC<4, NperDim>::weights;
    template<typename F, typename... Args>
    static double Integrate(const F& f, Args... args) {
        double sum = 0.0;
        #pragma omp parallel for reduction(+:sum)

        for (size_t i1 = 0; i1 < NperDim; i1++) {
            double x1 = 0.5 + 0.5 * points[i1];
            double w1 = weights[i1] * 0.5;

            for (size_t i2 = 0; i2 < NperDim; i2++) {
                double x2 = 0.5 + 0.5 * points[i2];
                double w2 = weights[i2] * 0.5;

                for (size_t i3 = 0; i3 < NperDim; i3++) {
                    double x3 = 0.5 + 0.5 * points[i3];
                    double w3 = weights[i3] * 0.5;

                    for (size_t i4 = 0; i4 < NperDim; i4++) {
                        double x4 = 0.5 + 0.5 * points[i4];
                        double w4 = weights[i4] * 0.5;
                        std::array<double, 4> x = {x1, x2, x3, x4};
                        sum += w1 * w2 * w3 * w4 * f(x, args...);
                    }
                }
            }
        }

        return sum;
    }
};

// SPECIALIZATION for 5D
template<unsigned NperDim>
struct QUAD<5, NperDim> : QUAD_STATIC<5, NperDim> {
    using QUAD_STATIC<5, NperDim>::points;
    using QUAD_STATIC<5, NperDim>::weights;
    template<typename F, typename... Args>
    static double Integrate(const F& f, Args... args) {
        double sum = 0.0;
        #pragma omp parallel for reduction(+:sum) collapse(3)

        for (size_t i1 = 0; i1 < NperDim; i1++) {
            double x1 = 0.5 + 0.5 * points[i1];
            double w1 = weights[i1] * 0.5;

            for (size_t i2 = 0; i2 < NperDim; i2++) {
                double x2 = 0.5 + 0.5 * points[i2];
                double w2 = weights[i2] * 0.5;

                for (size_t i3 = 0; i3 < NperDim; i3++) {
                    double x3 = 0.5 + 0.5 * points[i3];
                    double w3 = weights[i3] * 0.5;

                    for (size_t i4 = 0; i4 < NperDim; i4++) {
                        double x4 = 0.5 + 0.5 * points[i4];
                        double w4 = weights[i4] * 0.5;

                        for (size_t i5 = 0; i5 < NperDim; i5++) {
                            double x5 = 0.5 + 0.5 * points[i5];
                            double w5 = weights[i5] * 0.5;
                            std::array<double, 5> x = {x1, x2, x3, x4, x5};
                            sum += w1 * w2 * w3 * w4 * w5 * f(x, args...);
                        }
                    }
                }
            }
        }

        return sum;
    }
};

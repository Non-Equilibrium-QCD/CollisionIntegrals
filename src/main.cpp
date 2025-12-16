#include "Integral.cpp"
#include "IntegralRef.cpp"
#include <cstring>

// Define which implementation to use
// 0 = Original (Integral.cpp)
// 1 = Reference-style (IntegralRef.cpp)
#ifndef USE_REF_IMPL
#define USE_REF_IMPL 0
#endif

int main (int argc, char *argv[]) {
    // Check command line for --ref flag
    bool useRef = USE_REF_IMPL;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--ref") == 0) {
            useRef = true;
        } else if (std::strcmp(argv[i], "--orig") == 0) {
            useRef = false;
        }
    }

    if (useRef) {
        fmt::println(stderr, "Using REFERENCE implementation (IntegralRef.cpp)");
        IntegrateQCDRef::Setup();
        auto integrand = CollisionIntegralQCDRef::CollisionIntegral<gg_to_gg>;
        IntegrateQCDRef::Compute<gg_to_gg>(integrand);
    } else {
        fmt::println(stderr, "Using ORIGINAL implementation (Integral.cpp)");
        IntegrateQCD::Setup();
        auto integrand = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>;
        IntegrateQCD::Compute<gg_to_gg>(integrand);
    }

    return 0;
}

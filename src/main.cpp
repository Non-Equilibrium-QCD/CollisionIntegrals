#include "Integral.cpp"

int main (int argc, char *argv[]) {
    IntegrateQCD::Setup();

    // Instantiate collision integral for gg -> gg process
    auto integrand = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>;
    IntegrateQCD::Compute<gg_to_gg>(integrand);
    return 0;
}

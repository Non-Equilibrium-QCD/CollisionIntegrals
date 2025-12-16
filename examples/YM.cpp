#include "../src/Integral.cpp"

int main(int argc, char *argv[]) {

    fmt::println(stderr, "Using ORIGINAL implementation (Integral.cpp)");
    IntegrateQCD::Setup();
    auto integrand = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>;
    IntegrateQCD::Compute<gg_to_gg>(integrand);

    return 0;
}

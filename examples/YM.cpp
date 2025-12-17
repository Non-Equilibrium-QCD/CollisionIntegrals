#include "../src/Integral.cpp"
#include <fmt/color.h>

int main(int argc, char *argv[]) {
    (void)argc; (void)argv;

    std::string OutputFile = "OUTPUT/YM.dat";

    fmt::println(stderr, "Using ORIGINAL implementation (Integral.cpp)");
    IntegrateQCD::Setup();
    auto integrand = CollisionIntegralQCD::CollisionIntegral<gg_to_gg>;
    fmt::println(stderr, "Starting computation for gg->gg process...");
    fmt::println(stderr, "Output will be saved to {}",
                 fmt::styled(OutputFile, fmt::emphasis::bold | fmt::fg(fmt::color::green)));
    IntegrateQCD::Compute<gg_to_gg>(integrand, OutputFile);

    return 0;
}

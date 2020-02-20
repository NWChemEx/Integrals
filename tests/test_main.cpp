#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#include <tiledarray.h>

int main(int argc, char* argv[]) {
    // Initialize Everything and set world pointer
    auto& world = TA::initialize(argc, argv);

    // Run tests
    int res = Catch::Session().run(argc, argv);

    // Finalize Everything
    TA::finalize();

    return res;
}

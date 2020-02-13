#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#include <tiledarray.h>

// Pointer for world so it can be accessed by other test
TA::World* pworld;

int main(int argc, char* argv[]) {
    // Initialize Everything and set world pointer
    auto& world = TA::initialize(argc, argv);
    pworld = &world;

    // Run tests
    int res = Catch::Session().run(argc, argv);

    // Finalize Everything
    TA::finalize();

    return res;
}

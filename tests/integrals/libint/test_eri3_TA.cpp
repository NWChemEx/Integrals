#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <testing/testing.hpp>

TEST_CASE("ERI3C") {
    using integral_type = integrals::pt::eri3c<double>;
    const auto key      = "ERI3";

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = "h2o";
    const auto bs   = "sto-3g";
    auto mol        = testing::get_molecules().at(name);
    auto aos        = testing::get_bases().at(name).at(bs);
    auto tensors    = testing::get_data(world).at(name).at(bs);
    auto [X]        = mm.run_as<integral_type>(key, aos, aos, aos);
    REQUIRE(libchemist::ta_helpers::allclose(X, tensors.at("ERIs 3C")));
}

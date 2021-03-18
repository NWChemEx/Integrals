#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <testing/testing.hpp>

TEST_CASE("Kinetic") {
    using integral_type = integrals::pt::kinetic<double>;
    using size_vector   = integrals::type::size_vector;

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto name = "h2o";
    auto bs   = "sto-3g";
    auto aos  = testing::get_bases().at(name).at(bs);
    auto corr = testing::get_data(world).at(name).at(bs);
    mm.at("Kinetic").change_input("Tile size", size_vector{1, 2});
    auto [T]    = mm.at("Kinetic").run_as<integral_type>(aos, aos);
    auto T_corr = TA::retile(corr.at("Kinetic"), T.trange());
    REQUIRE(libchemist::ta_helpers::allclose(T, T_corr));
}

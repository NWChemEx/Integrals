#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <testing/testing.hpp>

TEST_CASE("ERI2C") {
    using integral_type = integrals::pt::eri2c<double>;
    using size_vector   = integrals::type::size_vector;

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto name = "h2o";
    auto bs   = "sto-3g";
    auto mol  = testing::get_molecules().at(name);
    auto aos  = testing::get_bases().at(name).at(bs);
    auto corr = testing::get_data(world).at(name).at(bs);

    mm.at("ERI2").change_input("Tile size", size_vector{1});
    auto [X]    = mm.at("ERI2").run_as<integral_type>(aos, aos);
    auto X_corr = TA::retile(corr.at("ERIs 2C"), X.trange());

    REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));
}

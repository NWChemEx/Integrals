#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <testing/testing.hpp>

TEST_CASE("STG2C") {
    using integral_type = integrals::pt::stg2c<double>;
    using size_vector   = integrals::type::size_vector;

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto name         = "h2o";
    auto bs           = "sto-3g";
    auto aos          = testing::get_bases().at(name).at(bs);
    auto corr         = testing::get_data(world).at(name).at(bs);
    auto stg_exponent = 1.0;
    mm.at("STG2").change_input("Tile size", size_vector{3});
    auto [X]    = mm.at("STG2").run_as<integral_type>(stg_exponent, aos, aos);
    auto corr_R = TA::retile(corr.at("STG 2C"), X.trange());
    REQUIRE(libchemist::ta_helpers::allclose(X, corr_R));
}
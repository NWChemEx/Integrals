#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <testing/testing.hpp>

TEST_CASE("STG3C") {
    using integral_type = integrals::pt::stg3c<double>;
    using size_vector   = integrals::type::size_vector;
    const auto key      = "STG3";

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name   = "h2o";
    const auto bs     = "sto-3g";
    auto mol          = testing::get_molecules().at(name);
    auto aos          = testing::get_bases().at(name).at(bs);
    auto stg_exponent = 1.0;

    auto [X]    = mm.run_as<integral_type>(key, stg_exponent, aos, aos, aos);
    auto corr_R = testing::get_data(world).at(name).at(bs).at("STG 3C");
    REQUIRE(libchemist::ta_helpers::allclose(X, corr_R));
}
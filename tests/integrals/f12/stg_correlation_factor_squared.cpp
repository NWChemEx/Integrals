#include "integrals/integralsmm.hpp"
#include "integrals/property_types.hpp"
#include <catch2/catch.hpp>
#include <libchemist/libchemist.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <libint2.hpp>
#include <testing/testing.hpp>

TEST_CASE("STG 4 Center Correlation Factor Squared") {
    using integral_type = integrals::pt::correlation_factor_squared_4c<double>;
    const auto key      = "STG 4 Center Correlation Factor Squared";

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto mols  = testing::get_molecules();
    auto bases = testing::get_bases();
    auto data  = testing::get_data(world);

    SECTION("H2") {
        const auto name = "h2";
        auto mol        = mols.at(name);
        for(auto bs : {"sto-3g", "cc-pvdz"}) {
            SECTION(bs) {
                auto aos     = bases.at(name).at(bs);
                auto tensors = data.at(name).at(bs);
                auto X_corr  = tensors.at("STG 4C correlation factor squared");
                auto [X] = mm.at(key).run_as<integral_type>(aos, aos, aos, aos);
                REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));
            }
        }
    }
}
#include "integrals/integralsmm.hpp"
#include "integrals/property_types.hpp"
#include <catch2/catch.hpp>
#include <libchemist/libchemist.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <libint2.hpp>
#include <testing/testing.hpp>

TEST_CASE("STG 2 Center Correlation Factor") {
    using integral_type = integrals::pt::correlation_factor_2c<double>;
    using overlap_type  = integrals::pt::overlap<double>;
    const auto key      = "STG 2 Center Correlation Factor";

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
                auto X_corr  = tensors.at("STG 2C correlation factor");
                auto [X]     = mm.at(key).run_as<integral_type>(aos, aos);
                REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));
            }
        }
    }
}

TEST_CASE("STG 4 Center Correlation Factor") {
    using integral_type = integrals::pt::correlation_factor_4c<double>;
    using overlap_type  = integrals::pt::overlap<double>;
    const auto key      = "STG 4 Center Correlation Factor";

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
                auto X_corr  = tensors.at("STG 4C correlation factor");
                auto [X]     = mm.at(key).run_as<integral_type>(aos, aos);
                REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));
            }
        }
    }
    SECTION("H2O") {
        const auto name = "h2o";
        auto mol        = mols.at(name);
        const auto bs   = "sto-3g";
        auto aos        = bases.at(name).at(bs);
        auto tensors    = data.at(name).at(bs);
        auto X_corr     = tensors.at("STG 4C correlation factor");

        auto [X] = mm.at(key).run_as<integral_type>(aos, aos, aos, aos);
        // REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));
    }
}
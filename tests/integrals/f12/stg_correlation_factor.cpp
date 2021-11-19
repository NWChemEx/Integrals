#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <chemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>

using namespace mokup;

TEST_CASE("STG 2 Center Correlation Factor") {
    using op_type       = simde::type::el_el_stg;
    using integral_type = simde::AOTensorRepresentation<2, op_type>;
    const auto key      = "STG2";

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2;
    const auto prop = property::stg_correlation_factor;

    chemist::Electron e;
    op_type stg;

    for(const auto& bs : {basis_set::sto3g, basis_set::ccpvdz}) {
        std::vector<mokup::basis_set> bs_key(2, bs);
        SECTION(as_string(name, bs)) {
            auto aos    = get_bases(name, bs);
            auto X_corr = get_ao_data(name, bs_key, prop, world);
            auto [X]    = mm.at(key).run_as<integral_type>(aos, stg, aos);
            REQUIRE(chemist::tensor::allclose(X, X_corr));
        }
    }
}

TEST_CASE("STG 4 Center Correlation Factor") {
    using op_type       = simde::type::el_el_stg;
    using integral_type = simde::AOTensorRepresentation<4, op_type>;
    const auto key      = "STG4";

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2;
    const auto prop = property::stg_correlation_factor;
    chemist::Electron e;
    op_type stg;

    for(const auto& bs : {basis_set::sto3g, basis_set::ccpvdz}) {
        std::vector<basis_set> bs_key(4, bs);
        SECTION(as_string(name, bs)) {
            auto aos    = get_bases(name, bs);
            auto X_corr = get_ao_data(name, bs_key, prop, world);
            auto [X] =
              mm.at(key).run_as<integral_type>(aos, aos, stg, aos, aos);
            REQUIRE(chemist::tensor::allclose(X, X_corr));
        }
    }
}

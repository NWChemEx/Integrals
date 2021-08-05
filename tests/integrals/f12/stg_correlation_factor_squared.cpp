#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>

TEST_CASE("STG 4 Center Correlation Factor Squared") {
    using op_type       = simde::type::el_el_stg;
    using integral_type = simde::AOTensorRepresentation<4, op_type>;
    const auto key      = "STG4";

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = mokup::molecule::h2;
    libchemist::Electron e;
    libchemist::STG stg;
    op_type f12(stg * stg, e, e);

    for(const auto& bs : {mokup::basis_set::sto3g, mokup::basis_set::ccpvdz}) {
        std::vector<mokup::basis_set> bs_key(4, bs);
        SECTION(as_string(name) + "/" + as_string(bs)) {
            auto aos     = mokup::get_bases().at(name).at(bs);
            auto tensors = mokup::get_ao_data(world).at(name).at(bs_key);
            auto X_corr =
              tensors.at(mokup::property::stg_correlation_factor_squared);
            auto [X] =
              mm.at(key).run_as<integral_type>(aos, aos, f12, aos, aos);
            REQUIRE(libchemist::tensor::allclose(X, X_corr));
        }
    }
}

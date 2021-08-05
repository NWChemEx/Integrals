#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/libchemist.hpp>
#include <libchemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>
#include <simde/tensor_representation/tensor_representation.hpp>

TEST_CASE("STG 4 Center dfdr Squared") {
    using op_type       = simde::type::el_el_f12_commutator;
    using integral_type = simde::AOTensorRepresentation<4, op_type>;
    const auto key      = "STG 4 Center dfdr Squared";

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = mokup::molecule::h2;
    const auto bs   = mokup::basis_set::sto3g;

    op_type fTf;

    for(const auto& bs : {mokup::basis_set::sto3g, mokup::basis_set::ccpvdz}) {
        std::vector<mokup::basis_set> bs_key(4, bs);
        SECTION(as_string(name) + "/" + as_string(bs)) {
            auto aos     = mokup::get_bases().at(name).at(bs);
            auto tensors = mokup::get_ao_data(world).at(name).at(bs_key);
            auto X_corr  = tensors.at(mokup::property::stg_dfdr_squared);

            auto [X] =
              mm.at(key).run_as<integral_type>(aos, aos, fTf, aos, aos);
            REQUIRE(libchemist::tensor::allclose(X, X_corr));
        }
    }
}

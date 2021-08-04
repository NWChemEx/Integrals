#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>

TEST_CASE("ERI4C CS") {
    using integral_type = simde::ERI4;

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    mm.change_input("ERI4 CS", "Screening Threshold", 0.005);

    const auto name = mokup::molecule::h2o;
    const auto bs   = mokup::basis_set::sto3g;
    auto mol        = mokup::get_molecules().at(name);
    auto aos        = mokup::get_bases().at(name).at(bs);
    std::vector bases{bs, bs, bs, bs};
    auto tensors = mokup::get_ao_data(world).at(name).at(bases);
    auto corr_S  = tensors.at(mokup::property::screened_eris);
    simde::type::el_el_coulomb r12;
    auto [X] = mm.at("ERI4 CS").run_as<integral_type>(aos, aos, r12, aos, aos);

    REQUIRE(libchemist::tensor::allclose(X, corr_S));
}

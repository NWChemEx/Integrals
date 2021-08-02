#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>

TEST_CASE("Nuclear") {
    using op_type       = simde::type::el_nuc_coulomb;
    using integral_type = simde::AOTensorRepresentation<2, op_type>;

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    auto name = mokup::molecule::h2o;
    auto bs   = mokup::basis_set::sto3g;
    auto mol  = mokup::get_molecules().at(name);
    auto aos  = mokup::get_bases().at(name).at(bs);
    std::vector bases{bs, bs};
    auto corr = mokup::get_ao_data(world).at(name).at(bases);

    // mm.at("Nuclear").change_input("Tile size", size_vector{6, 1});
    op_type riA(libchemist::Electron{}, mol);
    auto [V] = mm.at("Nuclear").run_as<integral_type>(aos, riA, aos);
    REQUIRE(libchemist::tensor::allclose(V, corr.at(mokup::property::nuclear)));
}

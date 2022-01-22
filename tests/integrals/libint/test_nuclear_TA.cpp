#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>
#include <tensorwrapper/tensor/allclose.hpp>

using namespace mokup;

TEST_CASE("Nuclear") {
    using op_type       = simde::type::el_nuc_coulomb;
    using integral_type = simde::AOTensorRepresentation<2, op_type>;

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    auto name = molecule::h2o;
    auto bs   = basis_set::sto3g;
    auto mol  = get_molecule(name);
    auto aos  = get_bases(name, bs);
    std::vector bases{bs, bs};
    auto corr = get_ao_data(name, bases, property::nuclear, world);

    // mm.at("Nuclear").change_input("Tile size", size_vector{6, 1});
    op_type riA(chemist::Electron{}, mol);
    auto [V] = mm.at("Nuclear").run_as<integral_type>(aos, riA, aos);
    REQUIRE(tensorwrapper::tensor::allclose(V, corr));
}

#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>
#include <tensorwrapper/tensor/allclose.hpp>

using namespace mokup;

TEST_CASE("ERI2C") {
    using op            = simde::type::el_el_coulomb;
    using integral_type = simde::AOTensorRepresentation<2, op>;
    using size_vector   = std::vector<std::size_t>;

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    auto name = molecule::h2o;
    auto bs   = basis_set::sto3g;
    auto mol  = get_molecule(name);
    auto aos  = get_bases(name, bs);
    std::vector bases{bs, bs};
    auto corr = get_ao_data(name, bases, property::eris, world);

    // mm.at("ERI2").change_input("Tile size", size_vector{1});
    simde::type::el_el_coulomb r12;
    auto [X] = mm.at("ERI2").run_as<integral_type>(aos, r12, aos);
    REQUIRE(tensorwrapper::tensor::allclose(X, corr));
}

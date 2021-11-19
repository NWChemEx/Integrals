#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <chemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>

using namespace mokup;

TEST_CASE("Yukawa2C") {
    using op_type       = simde::type::el_el_yukawa;
    using integral_type = simde::AOTensorRepresentation<2, op_type>;

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    auto name = molecule::h2o;
    auto bs   = basis_set::sto3g;
    auto aos  = get_bases(name, bs);
    std::vector bases{bs, bs};
    auto corr = get_ao_data(name, bases, property::yukawa, world);

    chemist::Electron e;
    op_type gr(chemist::operators::STG(1.0, 1.0), e, e);

    auto [X] = mm.at("Yukawa2").run_as<integral_type>(aos, gr, aos);
    REQUIRE(chemist::tensor::allclose(X, corr));
}

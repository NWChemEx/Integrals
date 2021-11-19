#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <chemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>

using namespace mokup;

TEST_CASE("STG3C") {
    using op_type       = simde::type::el_el_stg;
    using integral_type = simde::AOTensorRepresentation<3, op_type>;

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    auto name = molecule::h2o;
    auto bs   = basis_set::sto3g;
    auto aos  = get_bases(name, bs);
    std::vector bases{bs, bs, bs};
    auto corr = get_ao_data(name, bases, property::stg, world);
    chemist::Electron e;
    op_type stg(chemist::operators::STG(1.0, 1.0));
    auto [X] = mm.at("STG3").run_as<integral_type>(aos, stg, aos, aos);
    REQUIRE(chemist::tensor::allclose(X, corr));
}

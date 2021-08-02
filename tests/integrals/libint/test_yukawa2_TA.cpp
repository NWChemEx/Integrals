#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>

TEST_CASE("Yukawa2C") {
    using op_type       = simde::type::el_el_yukawa;
    using integral_type = simde::AOTensorRepresentation<2, op_type>;

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    auto name = mokup::molecule::h2o;
    auto bs   = mokup::basis_set::sto3g;
    auto aos  = mokup::get_bases().at(name).at(bs);
    std::vector bases{bs, bs};
    auto corr = mokup::get_ao_data(world).at(name).at(bases);

    libchemist::Electron e;
    op_type gr(libchemist::STG(1.0, 1.0), e, e);

    auto [X] = mm.at("Yukawa2").run_as<integral_type>(aos, gr, aos);
    REQUIRE(libchemist::tensor::allclose(X, corr.at(mokup::property::yukawa)));
}

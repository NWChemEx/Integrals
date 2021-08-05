#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>

TEST_CASE("ERI4C") {
    using op_type       = simde::type::el_el_coulomb;
    using integral_type = simde::AOTensorRepresentation<4, op_type>;

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = mokup::molecule::h2o;
    const auto bs   = mokup::basis_set::sto3g;
    auto mol        = mokup::get_molecules().at(name);
    auto aos        = mokup::get_bases().at(name).at(bs);
    std::vector bases{bs, bs, bs, bs};
    auto tensors = mokup::get_ao_data(world).at(name).at(bases);

    op_type r12;
    auto [X] = mm.at("ERI4").run_as<integral_type>(aos, aos, r12, aos, aos);
    libchemist::type::tensor corr(tensors.at(mokup::property::eris));
    REQUIRE(libchemist::tensor::allclose(X, corr));
}

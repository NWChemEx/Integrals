#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <chemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>

using namespace mokup;

TEST_CASE("ERI3C") {
    using op_type       = simde::type::el_el_coulomb;
    using integral_type = simde::AOTensorRepresentation<3, op_type>;
    const auto key1     = "ERI3";
    const auto key2     = "ERI3 CS";

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto mol        = get_molecule(name);
    auto aos        = get_bases(name, bs);
    std::vector bases{bs, bs, bs};
    auto corr = get_ao_data(name, bases, property::eris, world);
    op_type r12;
    auto [X] = mm.run_as<integral_type>(key1, aos, r12, aos, aos);
    REQUIRE(chemist::tensor::allclose(X, corr));
}

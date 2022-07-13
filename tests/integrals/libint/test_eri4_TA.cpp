#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>
#include <tensorwrapper/tensor/allclose.hpp>

using namespace mokup;

TEST_CASE("ERI4C") {
    using op_type       = simde::type::el_el_coulomb;
    using integral_type = simde::AOTensorRepresentation<4, op_type>;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto mol        = get_molecule(name);
    auto aos        = get_bases(name, bs);
    std::vector bases{bs, bs, bs, bs};
    auto corr = get_ao_data(name, bases, property::eris);

    op_type r12;
    auto [X] = mm.at("ERI4").run_as<integral_type>(aos, aos, r12, aos, aos);
    REQUIRE(tensorwrapper::tensor::allclose(X, corr));
}

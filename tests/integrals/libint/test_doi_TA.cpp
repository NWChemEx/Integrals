#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <chemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>

using namespace mokup;

TEST_CASE("DOI") {
    using op_type       = simde::type::el_el_delta;
    using integral_type = simde::EDOI;

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto mol        = get_molecule(name);
    auto aos        = get_bases(name, bs);
    std::vector bases{bs, bs};
    auto corr = get_ao_data(name, bases, property::dois, world);
    op_type d;
    auto [X] = mm.at("DOI").run_as<integral_type>(aos, d, aos);
    REQUIRE(chemist::tensor::allclose(X, corr));
}

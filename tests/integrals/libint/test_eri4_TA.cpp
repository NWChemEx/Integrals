#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/tensors/allclose.hpp>
#include <testing/testing.hpp>

using namespace testing;

TEST_CASE("ERI4C") {
    using integral_type = integrals::pt::eri4c<double>;

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto mol        = testing::get_molecules().at(name);
    auto aos        = testing::get_bases().at(name).at(bs);
    std::vector bases{bs, bs, bs, bs};
    auto tensors = testing::get_ao_data(world).at(name).at(bases);

    simde::type::el_el_coulomb r12;
    auto [X] = mm.at("ERI4").run_as<integral_type>(aos, aos, r12, aos, aos);
    libchemist::tensor corr(tensors.at(property::eris)));
    REQUIRE(libchemist::tensors::allclose(X, tcorr));
}

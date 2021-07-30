#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/tensors/allclose.hpp>
#include <testing/testing.hpp>

using namespace testing;

TEST_CASE("ERI3C") {
    using integral_type = integrals::pt::eri3c<double>;
    const auto key1     = "ERI3";
    const auto key2     = "ERI3 CS";

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto mol        = testing::get_molecules().at(name);
    auto aos        = testing::get_bases().at(name).at(bs);
    std::vector bases{bs, bs, bs};
    auto tensors = testing::get_ao_data(world).at(name).at(bases);
    simde::type::el_el_coulomb r12;
    auto [X] = mm.run_as<integral_type>(key1, aos, r12, aos, aos);
    libchemist::tensor corr(tensors.at(property::eris));
    REQUIRE(libchemist::ta_helpers::allclose(X, corr));
}

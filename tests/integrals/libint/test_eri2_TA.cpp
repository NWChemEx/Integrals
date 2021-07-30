#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/tensors/allclose.hpp>
#include <testing/testing.hpp>

using namespace testing;

TEST_CASE("ERI2C") {
    using integral_type = integrals::pt::eri2c<double>;
    using size_vector   = integrals::type::size_vector;

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto name = molecule::h2o;
    auto bs   = basis_set::sto3g;
    auto mol  = testing::get_molecules().at(name);
    auto aos  = testing::get_bases().at(name).at(bs);
    std::vector bases{bs, bs};
    auto corr = testing::get_ao_data(world).at(name).at(bases);

    mm.at("ERI2").change_input("Tile size", size_vector{1});
    simde::type::el_el_coulomb r12;
    auto [X]    = mm.at("ERI2").run_as<integral_type>(aos, r12, aos);
    auto X_corr = TA::retile(corr.at(property::eris), X.trange());
    libchemist::tensor corr(X_corr);
    REQUIRE(libchemist::tensor::allclose(X, corr));
}

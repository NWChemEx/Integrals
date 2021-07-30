#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/tensors/allclose.hpp>
#include <testing/testing.hpp>

using namespace testing;

TEST_CASE("Kinetic") {
    using integral_type = integrals::pt::kinetic<double>;
    using size_vector   = integrals::type::size_vector;

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto name = molecule::h2o;
    auto bs   = basis_set::sto3g;
    auto aos  = testing::get_bases().at(name).at(bs);
    std::vector bases{bs, bs};
    auto corr = testing::get_ao_data(world).at(name).at(bases);
    mm.at("Kinetic").change_input("Tile size", size_vector{1, 2});
    simde::type::el_kinetic t;
    auto [T]    = mm.at("Kinetic").run_as<integral_type>(aos, t, aos);
    auto T_corr = TA::retile(corr.at(property::kinetic), T.trange());
    libchemist::tensor corr(T_corr);
    REQUIRE(libchemist::tensor::allclose(T, corr));
}

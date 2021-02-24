#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <testing/testing.hpp>

TEST_CASE("Overlap") {
    using integral_type = integrals::pt::overlap<double>;
    using size_vector   = integrals::type::size_vector;
    using pair_vector   = integrals::type::pair_vector;

    // TODO: Better names which describe what's being tested
    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = "h2o";
    const auto bs   = "sto-3g";
    auto mol        = testing::get_molecules().at(name);
    auto aos        = testing::get_bases().at(name).at(bs);
    auto tensors    = testing::get_data(world).at(name).at(bs);
    auto corr_S     = tensors.at("Overlap");

    SECTION("Run 1") {
        mm.at("Overlap").change_input("Tile size", size_vector{3, 1, 1});
        auto [S] = mm.at("Overlap").run_as<integral_type>(aos, aos);
        auto X   = TA::retile(corr_S, S.trange());
        REQUIRE(libchemist::ta_helpers::allclose(S, X));
    }

    SECTION("Run 2") {
        pair_vector atom_groups{{0, 1}, {1, 2}, {2, 3}};
        mm.at("Overlap").change_input("Tile size", size_vector{100});
        mm.at("Overlap").change_input("Atom Tile Groups", atom_groups);
        auto [S] = mm.at("Overlap").run_as<integral_type>(aos, aos);
        auto X   = TA::retile(corr_S, S.trange());
        REQUIRE(libchemist::ta_helpers::allclose(S, X));
    }
}

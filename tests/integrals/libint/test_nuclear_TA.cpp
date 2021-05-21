#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <testing/testing.hpp>

using namespace testing;

TEST_CASE("Nuclear") {
    using integral_type = integrals::pt::nuclear<double>;
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
    mm.at("Nuclear").change_input("Tile size", size_vector{6, 1});
    auto [V]    = mm.at("Nuclear").run_as<integral_type>(mol, aos, aos);
    auto corr_V = TA::retile(corr.at(property::nuclear), V.trange());
    REQUIRE(libchemist::ta_helpers::allclose(V, corr_V));
}

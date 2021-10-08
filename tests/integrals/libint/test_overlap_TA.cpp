#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>

using namespace mokup;

TEST_CASE("Overlap") {
    using op_type       = simde::type::el_identity;
    using integral_type = simde::AOTensorRepresentation<2, op_type>;

    // TODO: Better names which describe what's being tested
    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto mol        = get_molecule(name);
    auto aos        = get_bases(name, bs);
    std::vector bases{bs, bs};
    auto corr_S = get_ao_data(name, bases, property::overlap, world);

    SECTION("Run 1") {
        // mm.at("Overlap").change_input("Tile size", size_vector{3, 1, 1});
        op_type I;
        auto [S] = mm.at("Overlap").run_as<integral_type>(aos, I, aos);
        REQUIRE(libchemist::tensor::allclose(S, corr_S));
    }

    // SECTION("Run 2") {
    //     pair_vector atom_groups{{0, 1}, {1, 2}, {2, 3}};
    //     mm.at("Overlap").change_input("Tile size", size_vector{100});
    //     mm.at("Overlap").change_input("Atom Tile Groups", atom_groups);
    //     auto [S] = mm.at("Overlap").run_as<integral_type>(aos, aos);
    //     auto X   = TA::retile(corr_S, S.trange());
    //     REQUIRE(libchemist::ta_helpers::allclose(S, X));
    // }
}

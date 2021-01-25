#include "../../tensors/H2O_STO3G_OVLP.hpp"
#include "../../test_common_TA.hpp"
#include "integrals/integrals.hpp"

TEST_CASE("Overlap") {
    using integral_type = integrals::pt::overlap<double>;
    using size_vector   = integrals::type::size_vector;
    using pair_vector   = integrals::type::pair_vector;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();

    // TODO: Better names which describe what's being tested

    SECTION("Run 1") {
        mm.at("Overlap").change_input("Tile size", size_vector{3, 1, 1});
        auto [S] = mm.at("Overlap").run_as<integral_type>(bs, bs);
        TensorType corr_S(S.world(), S.trange(), corr);
        REQUIRE(libchemist::ta_helpers::allclose(S, corr_S));
    }

    SECTION("Run 2") {
        pair_vector atom_groups{{0, 1}, {1, 2}, {2, 3}};
        mm.at("Overlap").change_input("Tile size", size_vector{100});
        mm.at("Overlap").change_input("Atom Tile Groups", atom_groups);
        auto [S] = mm.at("Overlap").run_as<integral_type>(bs, bs);
        TensorType corr_S(S.world(), S.trange(), corr);
        REQUIRE(libchemist::ta_helpers::allclose(S, corr_S));
    }
}

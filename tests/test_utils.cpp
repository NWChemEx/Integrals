#include "test_common_TA.hpp"
#include <integrals/nwx_libint/nwx_libint.hpp>
#include <integrals/nwx_TA/nwx_TA_utils.hpp>

TEST_CASE("Test Utilities") {
    auto [molecule, sto] = make_molecule("sto-3g");
    auto adz = libchemist::apply_basis("aug-cc-pvdz", molecule);

    auto atoms = std::vector<libint2::Atom>{
        {8, 0.000000000000000, -0.143222342980786, 0.000000000000000},
        {1, 1.638033502034240, 1.136556880358410, 0.000000000000000},
        {1, -1.638033502034240, 1.136556880358410, 0.000000000000000},
        };
    libint2::BasisSet obs("sto-3g", atoms);
    obs.set_pure(true);

    auto LI_sto = nwx_libint::make_basis(sto);
    REQUIRE(LI_sto == obs);

    auto LI_adz = nwx_libint::make_basis(adz);
    auto LI_bsets = nwx_libint::make_basis_sets({sto, adz});
    REQUIRE(LI_bsets == std::vector<libint2::BasisSet>{LI_sto, LI_adz});

    auto sets_max_prims = nwx_libint::sets_max_nprims(LI_bsets);
    REQUIRE(sets_max_prims == 9);

    auto sets_max_l = nwx_libint::sets_max_l(LI_bsets);
    REQUIRE(sets_max_l == 2);

    auto bs_range1 = nwx_TA::make_tiled_range(sto, 1);
    TA::TiledRange1 corr_range1{0, 5, 6, 7};
    REQUIRE(bs_range1 == corr_range1);

    auto bs_range2 = nwx_TA::make_tiled_range(sto, 6);
    TA::TiledRange1 corr_range2{0, 6, 7};
    REQUIRE(bs_range2 == corr_range2);

    auto bs_range3 = nwx_TA::make_tiled_range(sto, std::vector<std::size_t>{1});
    REQUIRE(bs_range3 == corr_range1);

    auto bs_range4 = nwx_TA::make_tiled_range(sto, std::vector<std::size_t>{5, 1, 1});
    REQUIRE(bs_range4 == corr_range1);

    auto bs_range5 = nwx_TA::make_tiled_range(sto, std::vector<std::size_t>{1, 1});
    REQUIRE(bs_range5 == corr_range1);

    auto bs_range6 = nwx_TA::make_tiled_range(sto, std::vector<std::size_t>{5, 1});
    REQUIRE(bs_range6 == corr_range1);

    auto i_range1 = nwx_TA::make_tiled_range(10, 2);
    TA::TiledRange1 corr_range3{0, 2, 4, 6, 8, 10};
    REQUIRE(i_range1 == corr_range3);

    auto i_range2 = nwx_TA::make_tiled_range(10, std::vector<std::size_t>{2});
    REQUIRE(i_range1 == corr_range3);

    auto i_range3 = nwx_TA::make_tiled_range(10, std::vector<std::size_t>{2, 2, 2, 2, 4});
    REQUIRE(i_range3 == corr_range3);

    auto shells = nwx_libint::aos2shells(LI_sto, 0, 5);
    std::vector<std::size_t> corr_shells{0, 1, 2};
    REQUIRE(shells == corr_shells);
    shells = nwx_TA::aos2shells(sto, 0, 5);
    REQUIRE(shells == corr_shells);

    auto trange1 = nwx_TA::make_trange({sto, sto}, std::vector<std::size_t>{1});
    std::vector<TA::TiledRange1> ranges{bs_range1, bs_range1};
    TA::TiledRange corr_trange1(ranges.begin(), ranges.end());
    REQUIRE(trange1 == corr_trange1);

    std::vector<std::vector<std::size_t>> atom_ranges = {{0, 1}, {2}};
    auto trange2 = nwx_TA::make_trange({sto, adz}, atom_ranges);
    TA::TiledRange1 adz_range{0, 32, 41};
    TA::TiledRange corr_trange2{corr_range2, adz_range};
    REQUIRE(trange2 == corr_trange2);

}
#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#include "test_common_TA.hpp"
#include "integrals/nwx_libint/nwx_libint.hpp"
#include "integrals/nwx_libint/nwx_libint_factory.hpp"
#include "integrals/nwx_TA/fill_2D_functor.hpp"
#include "integrals/nwx_TA/fill_3D_functor.hpp"
#include "integrals/nwx_TA/fill_4D_functor.hpp"

// Pointer for world so it can be accessed by other test
TA::World* pworld;

int main(int argc, char* argv[]) {
    // Initialize Everything and set world pointer
    auto& world = TA::initialize(argc, argv);
    pworld = &world;
    libint2::initialize();

    // Mock Params
    auto [molecule, bs] = make_molecule();
    auto tile_size = std::vector<std::size_t>{1};
    auto deriv = 0;
    auto thresh = 1.0E-16;
    auto origin = std::array<double, 3>{0,0,0};

    // Collect all basis sets
    std::vector<libchemist::AOBasisSet<double>> basis_sets{3, bs};

    // Details about multipole operator
//    auto nopers = libint2::operator_traits<libint2::Operator::emultipole3>::nopers;
//    libint2::BasisSet fake_basis{};
//    for (int i = 0; i < nopers; ++i) {
//        fake_basis.push_back(libint2::Shell::unit());
//    }
//    auto fake_range = nwx_TA::make_tiled_range(fake_basis, 1);

    // Make TA ranges, Libint basis sets, and Libint params based on basis sets
    std::vector<libint2::BasisSet> LIBasis_sets{};
    std::vector<TA::TiledRange1>   ranges{};
    std::size_t max_nprim = 0;
    int max_l = 0;

    for (auto i = 0; i < basis_sets.size(); ++i) {
        // Add tiled range based on each basis set
        ranges.push_back(nwx_TA::make_tiled_range(basis_sets[i], tile_size));

        // Make Libint basis set from LibChemist one
        LIBasis_sets.push_back(nwx_libint::make_basis(basis_sets[i]));

        // Find max_nprim and max_l over all basis sets
        auto max_nprim_i = libint2::max_nprim(LIBasis_sets[i]);
        auto max_l_i = libint2::max_l(LIBasis_sets[i]);
        max_nprim = std::max(max_nprim, max_nprim_i);
        max_l = std::max(max_l, max_l_i);
    }

    // Make the complete tiled range
    TA::TiledRange trange(ranges.begin(), ranges.end());

    // Make engine factory
    nwx_libint::LibintFactory<3, libint2::Operator::coulomb> factory(max_nprim, max_l, thresh, deriv);
    factory.mol = molecule;

    // Make TA fill functor
    nwx_TA::Fill3DFunctor<TA::TArrayD::value_type, libint2::Operator::coulomb> fill_a(LIBasis_sets, factory);
    nwx_TA::Fill3DFunctor<TiledArray::TSpArrayD::value_type, libint2::Operator::coulomb> fill_b(LIBasis_sets, factory);

    // Fill in the array, non-distributed
    double a_time_start = madness::wall_time();
    TA::TArrayD a(world, trange);
    for (auto it = begin(a); it != end(a); ++it) {
        auto tile = decltype(a)::value_type(a.trange().make_tile_range(it.index()));
        fill_a(tile, tile.range());
        *it = tile;
    }
    double a_time_stop = madness::wall_time();
    std::cout << a << std::endl;
    std::cout << a_time_stop - a_time_start << std::endl;

    // Fill in the array, distributed
    double b_time_start = madness::wall_time();
    auto b = TiledArray::make_array<TiledArray::TSpArrayD>(world, trange, fill_b);
    double b_time_stop = madness::wall_time();
    std::cout << b << std::endl;
    std::cout << b_time_stop - b_time_start << std::endl;

    // Potentially run other tests
    int res = Catch::Session().run(argc, argv);

    // Finalize Everything
    libint2::finalize();
    TA::finalize();

    return res;
}

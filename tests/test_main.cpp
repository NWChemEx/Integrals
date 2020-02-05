#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#include "test_common_TA.hpp"

// Pointer for world so it can be accessed by other test
TA::World* pworld;

int main(int argc, char* argv[]) {
    // Initialize Everything and set world pointer
    auto& world = TA::initialize(argc, argv);
    pworld = &world;
    libint2::initialize();

    // Mock Params
    auto bra1 = water_basis<double>();
    auto ket1 = water_basis<double>();
    auto bra2 = water_basis<double>();
    auto ket2 = water_basis<double>();
    auto tile_size = 1;
    auto deriv = 0;
    auto thresh = 1.0E-16;

    // Collect all basis sets
    std::vector<libchemist::AOBasisSet<double>> basis_sets{bra1, ket1, bra2, ket2};

    // Make TA ranges, Libint basis sets, and Libint params based on basis sets
    std::vector<libint2::BasisSet> LIBasis_sets{};
    std::vector<TA::TiledRange1>   ranges{};
    std::size_t max_nprim = 0;
    int max_l = 0;

    for (auto i = 0; i < basis_sets.size(); ++i) {
        // Add tiled range based on each basis set
        ranges.push_back(make_tiled_range(basis_sets[i], tile_size));

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

    // Non-distributed, better for testing/printing
    {
//        double a_time_start = madness::wall_time();
//        TA::TArrayD a(world, trange);
//        for (auto it = begin(a); it != end(a); ++it) {
//            auto tile = decltype(a)::value_type(a.trange().make_tile_range(it.index()));
//
//            // Make libint engine
//            // This is possibly questionable: doing this currently to get around issues with pointers
//            libint2::Engine l_engine(libint2::Operator::coulomb, max_nprim, max_l, thresh, deriv);
//            l_engine.set(libint2::BraKet::xx_xx);
//
//            // Grab bounds of current tile for shell mapping
//            auto bra1_lower = tile.range().lobound()[0];
//            auto ket1_lower = tile.range().lobound()[1];
//            auto bra1_upper = tile.range().upbound()[0];
//            auto ket1_upper = tile.range().upbound()[1];
//
//            auto bra2_lower = tile.range().lobound()[2];
//            auto ket2_lower = tile.range().lobound()[3];
//            auto bra2_upper = tile.range().upbound()[2];
//            auto ket2_upper = tile.range().upbound()[3];
//
//            // Initially set bra_offset to tile lower bound
//            auto bra1_offset = bra1_lower;
//
//            // Loop over Bra shells in tile
//            for (auto s1 : aos2shells(bra1, bra1_lower, bra1_upper)) {
//                // Number of basis functions in first shell
//                auto n1 = LIBasis_sets[0][s1].size();
//
//                // initially set ket_offset to tile lower bound
//                auto ket1_offset = ket1_lower;
//
//                // Loop over Ket shells in tile
//                for (auto s2 : aos2shells(ket1, ket1_lower, ket1_upper)) {
//                    // Number of basis functions in second shell
//                    auto n2 = LIBasis_sets[1][s2].size();
//
//                    auto bra2_offset = bra2_lower;
//
//                    for (auto s3 : aos2shells(bra2, bra2_lower, bra2_upper)) {
//                        auto n3 = LIBasis_sets[2][s3].size();
//
//                        auto ket2_offset = ket2_lower;
//
//                        for (auto s4 : aos2shells(ket2, ket2_lower, ket2_upper)) {
//                            auto n4 = LIBasis_sets[3][s4].size();
//
//                            // Compute integrals for current shells
//                            l_engine.compute(LIBasis_sets[0][s1], LIBasis_sets[1][s2],
//                                             LIBasis_sets[2][s3], LIBasis_sets[3][s4]);
//                            auto ints_shellset = l_engine.results()[0];
//
//                            // Default case of all zeros
//                            if (ints_shellset == nullptr) {
//                                std::fill(tile.begin(), tile.end(), 0.0);
//                                continue;
//                            }
//
//                            // Loop over basis functions in current shells
//                            for(auto f1 = 0ul; f1 != n1; ++f1) {
//                                for (auto f2 = 0ul; f2 != n2; ++f2) {
//                                    for (auto f3 = 0ul; f3 != n3; ++f3) {
//                                        for (auto f4 = 0ul; f4 != n4; ++f4) {
//                                            // Set index vector
//                                            std::vector<std::size_t> indexer{bra1_offset + f1, ket1_offset + f2,
//                                                                             bra2_offset + f3, ket2_offset + f4};
//
//                                            // Assign integral value
//                                            tile[indexer] = ints_shellset[(f1 * n2) + (f2 * n3) + (f3 * n4) + f4];
//                                        }
//                                    }
//                                }
//                            }
//                            ket2_offset += ket2.shell(s4).size();
//                        }
//                        bra2_offset += bra2.shell(s3).size();
//                    }
//                    // Increase offsets by finished shells' sizes
//                    ket1_offset += ket1.shell(s2).size();
//                }
//                bra1_offset += bra1.shell(s1).size();
//            }
//            *it = tile;
//        }
//        double a_time_stop = madness::wall_time();
//        std::cout << a << std::endl;
//        std::cout << a_time_stop - a_time_start << std::endl;
    }

    // Fill in the array, distributed
    double b_time_start = madness::wall_time();
    auto b = TiledArray::make_array<TiledArray::TSpArrayD>(world, trange,
        [&] (TiledArray::TSpArrayD::value_type& tile, const TiledArray::Range& range) -> float {
        // Initialize tile
        tile = TiledArray::TSpArrayD::value_type(range);

        // Make libint engine
        // This is possibly questionable: doing this currently to get around issues with pointers
        libint2::Engine tile_engine(libint2::Operator::coulomb, max_nprim, max_l, thresh, deriv);
        tile_engine.set(libint2::BraKet::xx_xx);

        // Grab bounds of current tile for shell mapping
        auto lowers = tile.range().lobound();
        auto uppers = tile.range().upbound();

        // Initially set the first offset to tile lower bound
        auto basis1_offset = lowers[0];
        for (auto s1 : aos2shells(LIBasis_sets[0], lowers[0], uppers[0])) {

            // Number of basis functions in first shell
            auto n1 = LIBasis_sets[0][s1].size();

            // initially set the second offset to tile lower bound
            auto basis2_offset = lowers[1];
            for (auto s2 : aos2shells(LIBasis_sets[1], lowers[1], uppers[1])) {

                // Number of basis functions in second shell
                auto n2 = LIBasis_sets[1][s2].size();

                // initially set the third offset to tile lower bound
                auto basis3_offset = lowers[2];
                for (auto s3 : aos2shells(LIBasis_sets[2], lowers[2], uppers[2])) {

                    // Number of basis functions in third shell
                    auto n3 = LIBasis_sets[2][s3].size();

                    // initially set the fourth offset to tile lower bound
                    auto basis4_offset = lowers[3];
                    for (auto s4 : aos2shells(LIBasis_sets[3], lowers[3], uppers[3])) {

                        // Number of basis functions in fourth shell
                        auto n4 = LIBasis_sets[3][s4].size();

                        // Compute integrals for current shells
                        tile_engine.compute(LIBasis_sets[0][s1],
                                            LIBasis_sets[1][s2],
                                            LIBasis_sets[2][s3],
                                            LIBasis_sets[3][s4]);
                        auto ints_shellset = tile_engine.results()[0];

                        // Default case of all zeros
//                        if (ints_shellset == nullptr) {
//                            std::fill(tile.begin(), tile.end(), 0.0);
//                            continue;
//                        }

                        // Loop over basis functions in current shells
                        for(auto f1 = 0ul; f1 != n1; ++f1) {
                            for (auto f2 = 0ul; f2 != n2; ++f2) {
                                for (auto f3 = 0ul; f3 != n3; ++f3) {
                                    for (auto f4 = 0ul; f4 != n4; ++f4) {
                                        // Set index vector
                                        std::vector<std::size_t> indexer{basis1_offset + f1, basis2_offset + f2,
                                                                         basis3_offset + f3, basis4_offset + f4};

                                        // Assign integral value
                                        tile[indexer] = ints_shellset[(f1 * n2) + (f2 * n3) + (f3 * n4) + f4];
                                    }
                                }
                            }
                        }
                        // Increase offsets by finished shells' sizes
                        basis4_offset += ket2.shell(s4).size();
                    }
                    basis3_offset += bra2.shell(s3).size();
                }
                basis2_offset += ket1.shell(s2).size();
            }
            basis1_offset += bra1.shell(s1).size();
        }
        // Return norm for new tile
        return tile.norm();
    });
    double b_time_stop = madness::wall_time();
    std::cout << b << std::endl;
    std::cout << b_time_stop - b_time_start << std::endl;

    /*
     * Stuff this function needs to have access to:
     * -Templated on appropriate value_type for tile
     * -Libint basis sets
     *      -Shells to which a range of AOs belong
     *      -Shell sizes
     *      -Shells themselves for the compute call
     *-Parameters to build the Libint engine
     *      -the operator with be specified by the module
     *      -the BraKet value?
     *      -max_nprim and max_l come from basis info
     *      -thresh and deriv are module inputs
     */

    // Potentially run other tests
    int res = Catch::Session().run(argc, argv);

    // Finalize Everything
    libint2::finalize();
    TA::finalize();

    return res;
}

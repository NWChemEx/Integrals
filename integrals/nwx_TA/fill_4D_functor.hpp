#pragma once
#include <libint2.hpp>
#include <tiledarray.h>
#include "integrals/nwx_libint/nwx_libint_factory.hpp"
#include "integrals/nwx_TA/nwx_TA_utils.hpp"

namespace nwx_TA {

    template<typename val_type, libint2::Operator op>
    struct Fill4DFunctor {
        using basis = libint2::BasisSet;

        std::vector<basis> LIBasis_sets;
        nwx_libint::LibintFactory<4, op> factory;

        Fill4DFunctor(std::vector<basis> LIBasis_sets, nwx_libint::LibintFactory<4, op> factory) :
                                  LIBasis_sets{std::move(LIBasis_sets)}, factory{std::move(factory)} {}

        float operator()(val_type& tile, const TiledArray::Range& range) {
            return _fill(tile, range);
        }

    private:
        float _fill(val_type& tile, const TiledArray::Range& range) {
            tile = val_type(range);

            // Make libint engine
            auto tile_engine = factory();

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

                            // Loop over basis functions in current shells
                            for(auto f1 = 0ul; f1 != n1; ++f1) {
                                for (auto f2 = 0ul; f2 != n2; ++f2) {
                                    for (auto f3 = 0ul; f3 != n3; ++f3) {
                                        for (auto f4 = 0ul; f4 != n4; ++f4) {
                                            // Set index vector
                                            std::vector<std::size_t> indexer{basis1_offset + f1, basis2_offset + f2,
                                                                             basis3_offset + f3, basis4_offset + f4};

                                            // Assign integral value
                                            if (ints_shellset == nullptr) {
                                                tile[indexer] = 0; // Default case of all zeroes
                                            } else {
                                                tile[indexer] = ints_shellset[(f1 * n2) + (f2 * n3) + (f3 * n4) + f4];
                                            }
                                        }
                                    }
                                }
                            }
                            // Increase offsets by finished shells' sizes
                            basis4_offset += LIBasis_sets[3][s4].size();
                        }
                        basis3_offset += LIBasis_sets[2][s3].size();
                    }
                    basis2_offset += LIBasis_sets[1][s2].size();
                }
                basis1_offset += LIBasis_sets[0][s1].size();
            }
            // Return norm for new tile
            return tile.norm();
        }
    };

} // namespace nwx_TA
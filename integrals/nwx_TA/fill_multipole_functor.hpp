#pragma once
#include <libint2.hpp>
#include <tiledarray.h>
#include "integrals/nwx_libint/nwx_libint_factory.hpp"
#include "integrals/nwx_TA/nwx_TA_utils.hpp"

namespace nwx_TA {

    template<typename val_type, libint2::Operator op>
    struct FillMultipoleFunctor {
        using basis = libint2::BasisSet;

        std::vector<basis> LIBasis_sets;
        nwx_libint::LibintFactory<2, op> factory;

        FillMultipoleFunctor(std::vector<basis> LIBasis_sets, nwx_libint::LibintFactory<2, op> factory) :
                LIBasis_sets{std::move(LIBasis_sets)}, factory{std::move(factory)} {}

        float operator()(val_type& tile, const TiledArray::Range& range) {
            return _fill(tile, range);
        }

    private:
        float _fill(val_type& tile, const TiledArray::Range& range) {
            tile = val_type(range);
            auto nopers = libint2::operator_traits<op>::nopers;

            // Make libint engine
            auto tile_engine = factory();
            const auto& buf_vec = tile_engine.results();

            // Grab bounds of current tile for shell mapping
            auto lowers = tile.range().lobound();
            auto uppers = tile.range().upbound();

            // Initially set the first offset to tile lower bound
            auto basis1_offset = lowers[1];
            for (auto s0 : aos2shells(LIBasis_sets[0], lowers[1], uppers[1])) {

                // Number of basis functions in first shell
                auto n1 = LIBasis_sets[0][s0].size();

                // initially set the second offset to tile lower bound
                auto basis2_offset = lowers[2];
                for (auto s1 : aos2shells(LIBasis_sets[1], lowers[2], uppers[2])) {

                    // Number of basis functions in second shell
                    auto n2 = LIBasis_sets[1][s1].size();

                    // Compute integrals for current shells
                    tile_engine.compute(LIBasis_sets[0][s0],
                                        LIBasis_sets[1][s1]);

                    // Loop over basis functions in current shells
                    for (auto c0 = 0ul; c0 != nopers; ++c0) {
                        auto ints_shellset = buf_vec[c0];
                        for (auto f1 = 0ul; f1 != n1; ++f1) {
                            for (auto f2 = 0ul; f2 != n2; ++f2) {
                                // Set index vector
                                std::vector<std::size_t> indexer{c0, basis1_offset + f1, basis2_offset + f2};

                                // Assign integral value
                                if (ints_shellset == nullptr) {
                                    tile[indexer] = 0; // Default case of all zeroes
                                } else {
                                    tile[indexer] = *(ints_shellset);
                                    ints_shellset++;
                                }
                            }
                        }
                    }
                    basis2_offset += LIBasis_sets[1][s1].size();
                }
                basis1_offset += LIBasis_sets[0][s0].size();
            }
            // Return norm for new tile
            return tile.norm();
        }
    };

} // namespace nwx_TA
#pragma once
#include <libint2.hpp>
#include <tiledarray.h>
#include "integrals/nwx_libint/nwx_libint_factory.hpp"
#include "integrals/nwx_TA/nwx_TA_utils.hpp"

namespace nwx_TA {

    template<typename val_type, libint2::Operator op, std::size_t NBases>
    struct FillNDFunctor {

        using basis = libint2::BasisSet;
        using size_vec = std::vector<std::size_t>;

        std::vector<basis> LIBasis_sets;
        nwx_libint::LibintFactory<NBases, op> factory;

        FillNDFunctor() = default;

        float operator()(val_type& tile, const TiledArray::Range& range) {
            return _fill(tile, range);
        }

    private:
        float _fill(val_type& tile, const TiledArray::Range& range) {
            tile = val_type(range);

            // Make libint engine
            auto tile_engine = factory();
            const auto& buf_vec = tile_engine.results();

            size_vec offsets(NBases);
            size_vec shells(NBases);

            _index_shells(tile, 0, shells, tile_engine, offsets, range);

            // Return norm for new tile
            return tile.norm();
        }

        void _index_shells(val_type& tile,
                           int depth,
                           size_vec& shells,
                           libint2::Engine& tile_engine,
                           size_vec& offsets,
                           const TiledArray::Range& range) {
            offsets[depth] = range.lobound()[depth];
            for (auto s : aos2shells(LIBasis_sets[depth], range.lobound()[depth], range.upbound()[depth])) {
                shells[depth] = s;

                // Compute integrals for current shells
                if (depth == (NBases - 1)) {
                    _call_libint(tile_engine, shells, std::make_index_sequence<NBases>());
                    auto ints_shellset = tile_engine.results()[0];

                    auto indexer = offsets;
                    int int_i = 0;
                    _fill_from_libint(tile, ints_shellset, 0, int_i, indexer, offsets, shells);
                } else {
                    _index_shells(tile, depth + 1, shells, tile_engine, offsets, range);
                }

                offsets[depth] += LIBasis_sets[depth][shells[depth]].size();
            }
        }

        void _fill_from_libint(val_type& tile,
                               const double* ints_shellset,
                               int depth,
                               int& int_i,
                               size_vec& indexer,
                               const size_vec& offsets,
                               const size_vec& shells) {
            for (auto f = 0ul; f != LIBasis_sets[depth][shells[depth]].size(); ++f) {
                // Assign integral value
                if (depth == (NBases - 1)) {
                    if (ints_shellset == nullptr) {
                        tile[indexer] = 0; // Default case of all zeroes
                    } else {
                        tile[indexer] = ints_shellset[int_i];
                        int_i++;
                    }
                } else {
                    _fill_from_libint(tile, ints_shellset, depth + 1, int_i,indexer, offsets, shells);
                }
                indexer[depth]++;
            }
            indexer[depth] = offsets[depth];
        }

        template<std::size_t... Is>
        void _call_libint(libint2::Engine& tile_engine, size_vec shells, std::index_sequence<Is ...>) {
            tile_engine.compute(LIBasis_sets[Is][shells[Is]] ...);
        }
    };

} // namespace nwx_TA
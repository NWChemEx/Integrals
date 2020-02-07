#pragma once
#include <libint2.hpp>
#include <tiledarray.h>
#include "integrals/types.hpp"
#include <libchemist/basis_set/ao_basis_set/ao_basis_set.hpp>

namespace nwx_TA {

    using size = integrals::type::size;

    // Make tile range based on basis set
    template<typename T = double>
    TA::TiledRange1 make_tiled_range(libchemist::AOBasisSet<T> basis_set, std::size_t tile_size) {
        std::vector<std::size_t> bounds{0};

        auto span = 0;
        for (auto ishell = 0; ishell < basis_set.n_shells(); ++ishell) {
            span += basis_set.shell(ishell).size();

            if (span < tile_size) continue;

            auto bound = bounds.back() + span;
            bounds.push_back(bound);
            span = 0;
        }
        if (span != 0) bounds.push_back(bounds.back() + span);

        return TA::TiledRange1(bounds.begin(), bounds.end());
    }

    template<typename T = double>
    TA::TiledRange1 make_tiled_range(libchemist::AOBasisSet<T> basis_set, std::vector<size> tile_sizes) {
        if (tile_sizes.size() == 1) {
            return make_tiled_range(basis_set, tile_sizes[0]);
        }

        std::vector<size> bounds{0};
        auto sizes_index = 0;

        auto span = 0;
        for (auto ishell = 0; ishell < basis_set.n_shells(); ++ishell) {
            span += basis_set.shell(ishell).size();

            if (span < tile_sizes[sizes_index]) continue;

            auto bound = bounds.back() + span;
            bounds.push_back(bound);
            span = 0;

            sizes_index++;
            if (sizes_index >= tile_sizes.size()) {
                sizes_index = 0;
            }
        }
        if (span != 0) bounds.push_back(bounds.back() + span);

        return TA::TiledRange1(bounds.begin(), bounds.end());
    }

    TA::TiledRange1 make_tiled_range(libint2::BasisSet basis_set, size tile_size) {
        std::vector<size> bounds{0};

        auto span = 0;
        for (auto ishell = 0; ishell < basis_set.size(); ++ishell) {
            span += basis_set[ishell].size();

            if (span < tile_size) continue;

            auto bound = bounds.back() + span;
            bounds.push_back(bound);
            span = 0;
        }
        if (span != 0) bounds.push_back(bounds.back() + span);

        return TA::TiledRange1(bounds.begin(), bounds.end());
    }

    TA::TiledRange1 make_tiled_range(libint2::BasisSet basis_set, std::vector<size> tile_sizes) {
        if (tile_sizes.size() == 1) {
            return make_tiled_range(basis_set, tile_sizes[0]);
        }

        std::vector<size> bounds{0};
        auto sizes_index = 0;

        auto span = 0;
        for (auto ishell = 0; ishell < basis_set.size(); ++ishell) {
            span += basis_set[ishell].size();

            if (span < tile_sizes[sizes_index]) continue;

            auto bound = bounds.back() + span;
            bounds.push_back(bound);
            span = 0;

            sizes_index++;
            if (sizes_index >= tile_sizes.size()) {
                sizes_index = 0;
            }
        }
        if (span != 0) bounds.push_back(bounds.back() + span);

        return TA::TiledRange1(bounds.begin(), bounds.end());
    }

    TA::TiledRange1 make_tiled_range(size upper, size tile_size) {
        std::vector<size> bounds{};

        for (auto bound = 0; bound < upper; bound+=tile_size) {
            bounds.push_back(bound);
        }
        if (bounds.back() != upper) bounds.push_back(upper);

        return TA::TiledRange1(bounds.begin(), bounds.end());
    }

    TA::TiledRange1 make_tiled_range(size upper, std::vector<size> tile_sizes) {
        if (tile_sizes.size() == 1) {
            return make_tiled_range(upper, tile_sizes[0]);
        }

        std::vector<size> bounds{};
        auto sizes_index = 0;
        auto bound = 0;

        while (bound < upper) {
            bounds.push_back(bound);

            bound += tile_sizes[sizes_index];

            sizes_index++;
            if (sizes_index >= tile_sizes.size()) {
                sizes_index = 0;
            }
        }
        if (bounds.back() != upper) bounds.push_back(upper);

        return TA::TiledRange1(bounds.begin(), bounds.end());
    }

    std::vector<size> aos2shells(libint2::BasisSet basis_set, size lower, size upper) {
        std::vector<size> return_vec;

        for (auto ishell = 0, offset = 0; ishell < basis_set.size(); ++ishell) {
            if (offset >= upper) break;
            if (offset >= lower) return_vec.push_back(ishell);
            offset += basis_set[ishell].size();
        }

        return return_vec;
    }

} // namespace nwx_TA
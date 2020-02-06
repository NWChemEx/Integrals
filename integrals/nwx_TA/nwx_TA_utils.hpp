#pragma once
#include <libint2.hpp>
#include <tiledarray.h>
#include <libchemist/basis_set/ao_basis_set/ao_basis_set.hpp>

namespace nwx_TA {

    // Make tile range based on BS
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

    std::vector<std::size_t> aos2shells(libint2::BasisSet basis_set, std::size_t lower, std::size_t upper) {
        std::vector<std::size_t> return_vec;

        for (auto ishell = 0, offset = 0; ishell < basis_set.size(); ++ishell) {
            if (offset >= upper) break;
            if (offset >= lower) return_vec.push_back(ishell);
            offset += basis_set[ishell].size();
        }

        return return_vec;
    }

} // namespace nwx_TA
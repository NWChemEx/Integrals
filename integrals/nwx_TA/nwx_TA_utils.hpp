#pragma once
#include <libint2.hpp>
#include <tiledarray.h>
#include "integrals/types.hpp"
#include <libchemist/basis_set/ao_basis_set/ao_basis_set.hpp>

namespace nwx_TA {

    using size = integrals::type::size;

    // Make tile range based on basis set
    TA::TiledRange1 make_tiled_range(libint2::BasisSet basis_set, size tile_size);
    TA::TiledRange1 make_tiled_range(libint2::BasisSet basis_set, std::vector<size> tile_sizes);
    TA::TiledRange1 make_tiled_range(size upper, size tile_size);
    TA::TiledRange1 make_tiled_range(size upper, std::vector<size> tile_sizes);

    std::vector<size> aos2shells(libint2::BasisSet basis_set, size lower, size upper);

} // namespace nwx_TA
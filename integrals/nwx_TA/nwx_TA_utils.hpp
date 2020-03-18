#pragma once
#include <libint2.hpp>
#include <tiledarray.h>
#include "integrals/types.hpp"
#include <libchemist/basis_set/ao_basis_set/ao_basis_set.hpp>

namespace nwx_TA {

    using size = integrals::type::size;

    /** @brief Make a TA TiledRange from basis set information
     *
     *  Takes a LibInt2 BasisSet @p basis_set and a tile size @p tile_size and
     *  returns a TiledArray TiledRange1 where the shells are tiled together
     *  sequentially and the tile sizes are greater than or equal to @p tile_size.
     *
     *  @param[in] bs The LibInt2 BasisSet to be tiled
     *  @param[in] tile_size The minimum size for the tiles
     *  @returns The TilesArray TiledRange1 for the basis set
     */
    TA::TiledRange1 make_tiled_range(libint2::BasisSet basis_set, size tile_size);

    /** @brief Make a TA TiledRange from basis set information
     *
     *  Takes a LibInt2 BasisSet @p basis_set and an std::vector of tile sizes
     *  @p tile_sizes, allowing for variable tile sizes. The tile sizes are looped
     *  over if end of the vector is reached before the entire basis set is tiled.
     *
     *  @param[in] bs The LibInt2 BasisSet to be tiled
     *  @param[in] tile_sizes The vector of minimum tiles sizes
     *  @returns The TilesArray TiledRange1 for the basis set
     */
    TA::TiledRange1 make_tiled_range(libint2::BasisSet basis_set, std::vector<size> tile_sizes);

    /** @brief Produces a TiledArray TiledRange1 from 0 to the specified upper value @p upper
     *         with tile size @p tile_size.
     *
     *  @param[in] upper The upper value of the range
     *  @param[in] tile_size The tile size of the tiles
     *  @returns The TilesArray TiledRange1 for the range
     */
    TA::TiledRange1 make_tiled_range(size upper, size tile_size);

    /** @brief Produces a TiledArray TiledRange1 from 0 to the specified upper value @p upper
     *         with variable tile sizes @p tile_sizes.
     *
     *  The tile sizes are looped over if the end of the vector is reached before the end of the range.
     *
     *  @param[in] upper The upper value of the range
     *  @param[in] tile_sizes The vector of tile sizes for the tiles
     *  @returns The TilesArray TiledRange1 for the range
     */
    TA::TiledRange1 make_tiled_range(size upper, std::vector<size> tile_sizes);

    /** @brief Find the shells that contain the specified AOs.
     *
     *  Given a LibInt2 BasisSet @p basis_set and a lower @p lower and upper @p upper
     *  AO index, returns a std::vector of shell indices that contain the AOs between
     *  @p lower and @p upper.
     *
     *  @param[in] basis_set The LibInt2 BasisSet containing the AOs
     *  @param[in] lower The lower value of the AO range
     *  @param[in] upper The upper value of the AO range
     *  @returns The std::vector of the shell indices
     */
    std::vector<size> aos2shells(libint2::BasisSet basis_set, size lower, size upper);

    /** @brief Make trange from basis set information.
     *
     *  Takes a std::vector of LibInt2 BasisSets @p basis_sets, a std::vector of
     *  tile sizes @p tile_sizes, and a std::vector of TiledArray TiledRange1 @p ranges
     *  and returns a TiledArray TiledRange composed of the any values in @p ranges
     *  and the tiled ranges of the basis sets.
     *
     *  @param[in] basis_sets The vector of LibInt2 basis sets
     *  @param[in] tile_sizes The vector of tile sizes
     *  @param[in] ranges A potentially empty vector of premade TiledRange1
     *  @returns The TiledArray TiledRange made from all of the TiledRange1
     */
    TA::TiledRange make_trange(const std::vector<libint2::BasisSet>& basis_sets,
                               const std::vector<size>& tile_sizes,
                               std::vector<TA::TiledRange1> ranges = {});

} // namespace nwx_TA
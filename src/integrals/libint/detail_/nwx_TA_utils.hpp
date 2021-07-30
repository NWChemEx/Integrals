#pragma once
#include "integrals/types.hpp"
#include <libchemist/basis_set/ao_basis_set.hpp>
#include <tiledarray.h>

namespace nwx_TA {

template<typename T>
using basis = libchemist::AOBasisSet<T>;
using size  = std::size_t;

/** @brief Make a TA TiledRange from basis set information
 *
 *  Takes a LibChemist BasisSet @p basis_set and a tile size @p tile_size and
 *  returns a TiledArray TiledRange1 where the shells are tiled together
 *  sequentially and the tile sizes are greater than or equal to @p tile_size.
 *
 *  @param[in] bs The LibChemist BasisSet to be tiled
 *  @param[in] tile_size The minimum size for the tiles
 *  @returns The TilesArray TiledRange1 for the basis set
 */
TA::TiledRange1 make_tiled_range(basis<double> basis_set, size tile_size);
TA::TiledRange1 make_tiled_range(basis<float> basis_set, size tile_size);

/** @brief Make a TA TiledRange from basis set information
 *
 *  Takes a LibChemist BasisSet @p basis_set and an std::vector of tile sizes
 *  @p tile_sizes, allowing for variable tile sizes. The tile sizes are looped
 *  over if end of the vector is reached before the entire basis set is tiled.
 *
 *  @param[in] bs The LibChemist BasisSet to be tiled
 *  @param[in] tile_sizes The vector of minimum tiles sizes
 *  @returns The TilesArray TiledRange1 for the basis set
 */
TA::TiledRange1 make_tiled_range(basis<double> basis_set,
                                 std::vector<size> tile_sizes);
TA::TiledRange1 make_tiled_range(basis<float> basis_set,
                                 std::vector<size> tile_sizes);

/** @brief Produces a TiledArray TiledRange1 from 0 to the specified upper value
 * @p upper with tile size @p tile_size.
 *
 *  @param[in] upper The upper value of the range
 *  @param[in] tile_size The tile size of the tiles
 *  @returns The TilesArray TiledRange1 for the range
 */
TA::TiledRange1 make_tiled_range(size upper, size tile_size);

/** @brief Produces a TiledArray TiledRange1 from 0 to the specified upper value
 * @p upper with variable tile sizes @p tile_sizes.
 *
 *  The tile sizes are looped over if the end of the vector is reached before
 * the end of the range.
 *
 *  @param[in] upper The upper value of the range
 *  @param[in] tile_sizes The vector of tile sizes for the tiles
 *  @returns The TilesArray TiledRange1 for the range
 */
TA::TiledRange1 make_tiled_range(size upper, std::vector<size> tile_sizes);

/** @brief Find the shells that contain the specified AOs.
 *
 *  Given a LibChemist BasisSet @p basis_set and a lower @p lower and upper @p
 * upper AO index, returns a std::vector of shell indices that contain the AOs
 * between
 *  @p lower and @p upper.
 *
 *  @param[in] basis_set The LibChemist BasisSet containing the AOs
 *  @param[in] lower The lower value of the AO range
 *  @param[in] upper The upper value of the AO range
 *  @returns The std::vector of the shell indices
 */
std::vector<size> aos2shells(basis<double> basis_set, size lower, size upper);
std::vector<size> aos2shells(basis<float> basis_set, size lower, size upper);

/** @brief Make trange from basis set information.
 *
 *  Takes a std::vector of LibChemist BasisSets @p basis_sets, a std::vector of
 *  tile sizes @p tile_sizes, and a std::vector of TiledArray TiledRange1 @p
 * ranges and returns a TiledArray TiledRange composed of the any values in @p
 * ranges and the tiled ranges of the basis sets.
 *
 *  @param[in] basis_sets The vector of LibChemist basis sets
 *  @param[in] tile_sizes The vector of tile sizes
 *  @param[in] ranges A potentially empty vector of premade TiledRange1
 *  @returns The TiledArray TiledRange made from all of the TiledRange1
 */
TA::TiledRange make_trange(const std::vector<basis<double>>& basis_sets,
                           const std::vector<size>& tile_sizes,
                           std::vector<TA::TiledRange1> ranges = {});

TA::TiledRange make_trange(const std::vector<basis<float>>& basis_sets,
                           const std::vector<size>& tile_sizes,
                           std::vector<TA::TiledRange1> ranges = {});

/** @brief Make trange from basis set information.
 *
 *  Takes a std::vector of LibChemist BasisSets @p basis_sets, a std::vector of
 *  atom ranges @p tile_sizes, and a std::vector of TiledArray TiledRange1 @p
 * ranges and returns a TiledArray TiledRange composed of any values in @p
 * ranges and the tiled ranges of the basis sets.
 *
 *  Allows for tiling by atom blocks
 *
 *  @param[in] basis_sets The vector of LibChemist basis sets
 *  @param[in] atom_ranges The vector of atoms ranges per tile
 *  @param[in] ranges A potentially empty vector of premade TiledRange1
 *  @returns The TiledArray TiledRange made from all of the TiledRange1
 */
TA::TiledRange make_trange(
  const std::vector<basis<double>>& basis_sets,
  const std::vector<std::pair<size, size>>& atom_ranges,
  std::vector<TA::TiledRange1> ranges = {});

TA::TiledRange make_trange(
  const std::vector<basis<float>>& basis_sets,
  const std::vector<std::pair<size, size>>& atom_ranges,
  std::vector<TA::TiledRange1> ranges = {});

/** @brief Make trange from basis set and tiling information.
 *
 *  Based on tiling settings, return the appropriate TiledRange.
 *
 *  @param[in] basis_sets The vector of LibChemist basis sets
 *  @param[in] tile_sizes The vector of tile sizes
 *  @param[in] atom_ranges The vector of atoms ranges per tile
 *  @param[in] ranges A potentially empty vector of premade TiledRange1
 *  @returns The TiledArray TiledRange made from all of the TiledRange1
 */
TA::TiledRange select_tiling(
  const std::vector<basis<double>>& basis_sets,
  const std::vector<size>& tile_sizes,
  const std::vector<std::pair<size, size>>& atom_ranges,
  std::vector<TA::TiledRange1> ranges = {});

TA::TiledRange select_tiling(
  const std::vector<basis<float>>& basis_sets,
  const std::vector<size>& tile_sizes,
  const std::vector<std::pair<size, size>>& atom_ranges,
  std::vector<TA::TiledRange1> ranges = {});

} // namespace nwx_TA

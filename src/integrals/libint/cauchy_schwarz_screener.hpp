#pragma once
#include "nwx_libint_factory.hpp"
#include <libint2.hpp>
#include <tiledarray.h>

namespace nwx_libint {

// Class to handle the Cauchy-Schwarz screening of an integral
template<std::size_t NBases>
struct CauchySchwarzScreener {
    using size_type  = integrals::type::size;
    using size_vec   = std::vector<size_type>;
    using basis_type = libint2::BasisSet;
    using basis_vec  = std::vector<basis_type>;
    using double_vec = std::vector<double>;
    using approx_vec = std::vector<double_vec>;

    // Matrices to hold the approximations for the two sides of the integral
    approx_vec cs_mat1 = {};
    approx_vec cs_mat2 = {};

    /** @brief Check if a whole tile passes screening
     *
     *  @param basis_sets The basis sets of the integral
     *  @param range The range of the tile
     *  @param cs_thresh The screening threshold
     *  @returns true if the tile is screened out
     */
    bool tile(const basis_vec& basis_sets, const TiledArray::Range& range,
              double cs_thresh);

    /** @brief Check if a shell set passes screening
     *
     *  @param shells The indices of the shells
     *  @param cs_thresh The screening threshold
     *  @returns true if the tile is screened out
     */
    bool shellset(size_vec shells, double cs_thresh);

    /** @brief Set input vectors as tile specific subset of this objects
     * screening vectors
     *
     *  @param basis_sets The basis sets of the integral
     *  @param range The range of the tile
     *  @param mat1 The vector that is a subset of cs_mat1
     *  @param mat2 The vector that is a subset of cs_mat2
     */
    void set_sub_screen(const basis_vec& basis_sets,
                        const TiledArray::Range& range, approx_vec& mat1,
                        approx_vec& mat2) const;

}; // Class CauchySchwarzScreener

extern template class CauchySchwarzScreener<2>;
extern template class CauchySchwarzScreener<3>;
extern template class CauchySchwarzScreener<4>;

} // namespace nwx_libint
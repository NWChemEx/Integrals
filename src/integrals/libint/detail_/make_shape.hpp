#pragma once
#include <libint2.hpp>
#include <simde/types.hpp>

namespace integrals::detail_ {

/** @brief Given a vector of basis sets, compute the shape of the corresponding
 *         integral.
 *
 *  @param[in] bases A vector of LibInt2 BasisSets
 *  @returns A unique_ptr for the resulting shape
 */
inline auto make_shape(const std::vector<libint2::BasisSet>& bases) {
    using shape_t   = typename simde::type::tensor::shape_type;
    using extents_t = typename shape_t::extents_type;

    /// Count up the extents from the shells in the basis sets
    extents_t extents{};
    for(auto& set : bases) {
        std::size_t set_extent = 0;
        for(auto& shell : set) { set_extent += shell.size(); }
        extents.push_back(set_extent);
    }

    /// Make the unique pointer to the shape
    return std::make_unique<shape_t>(extents);
}

} // namespace integrals::detail_
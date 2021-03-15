#pragma once
#include "integrals/types.hpp"
#include <property_types/ao_integrals/type_traits.hpp>
#include <property_types/ao_integrals/utilities/unpack_spaces.hpp>

namespace integrals {

/** @brief Given the inputs to a module unpacks the spaces for the integral.
 *
 *  This function wraps the logic required to unpack the spaces provided to an
 *  integral module. It can be used to unpack the spaces for normal AO integrals
 *  as well as for transformed variants.
 *
 *  @todo Should SpaceType be determined by template meta-programming? Similarly
 *        for the bra/ket prefixes?
 *
 *  @tparam PropType The property type of the integral. Assumed to be a valid
 *                   property type for an AO integral or a transformed variant
 *                   of one.
 *  @tparam SpaceType The type of the spaces being unpacked. For a normal AO
 *                    integral this is going to be libchemist::AOSpace, for
 *                    transformed variants its libchemist::BaseSpace.
 *
 *  @input[in] inputs The input map provided to the module. The spaces will be
 *                    unpacked from this object.
 *  @input[in] bra    The prefix of the key for bras. Default is `"bra"`.
 *  @input[in] ket    The prefix of the key for kets. Default is "ket"
 *
 *  @return An std::vector whose elements are spaces for the integrals. The
 *          orderis left to right across the braket, e.g., the return for the
 *          ERI integral (i|ab) would be the space for i, the space for a, and
 *          then the space for b.
 *
 *  @throw std::bad_alloc if there is insufficient memory to create the return.
 *                        Strong throw guarantee.
 *  @throw std::runtime_error if the inputs do not contain the correct keys.
 *                            Strong throw guarantee.
 */
template<typename PropType, typename SpaceType>
auto unpack_spaces(const sde::type::input_map& inputs,
                   const std::string& bra = "bra",
                   const std::string& ket = "ket") {
    using element_type = double; // TODO: Get from PropType via TMP
    using key_vector   = std::vector<std::string>;
    using namespace property_types::ao_integrals;

    constexpr auto n_centers = n_centers_v<PropType>;
    key_vector keys;

    const auto bra1 = bra + " 1", bra2 = bra + " 2";
    const auto ket1 = ket + " 1", ket2 = ket + " 2";

    if constexpr(is_doi_v<PropType>) {
        keys = key_vector{bra, bra, ket, ket};
    } else if constexpr(n_centers == 2) {
        keys = key_vector{bra, ket};
    } else if constexpr(n_centers == 3) {
        keys = key_vector{bra, ket1, ket2};
    } else if constexpr(n_centers == 4) {
        keys = key_vector{bra1, bra2, ket1, ket2};
    }

    std::vector<std::decay_t<SpaceType>> spaces;
    for(const auto& key_i : keys)
        spaces.push_back(inputs.at(key_i).value<SpaceType>());
    return spaces;
}

/** @brief Unpacks the AO basis sets from the inputs to an integral module.
 *
 *  This function is largely implemented by `unpack_spaces`. In addition to
 *  unpacking the spaces this function calls the `.basis_set()` member on each
 *  of the spaces and returns the result.
 *
 *  @tparam PropType The property type of the integral whose inputs are being
 *                   unpacked.
 *
 *  @param[in] inputs The input map provided to the module. The AO basis sets
 *                    will be unpacked from this object.
 *
 *  @return An std::vector whose elements are the basis sets for the integral.
 *          The order goes from left to right across the integral. For example
 *          the returns for the integral @f$(ij|ab)@f$ will be the AO basis sets
 *          for @f$i@f$, @f$j@f$, @f$a@f$, and then @f$b@f$.
 *
 *  @throw std::bad_alloc if there is insufficient memory to create the return.
 *                        Strong throw guarantee.
 *  @throw std::runtime_error if @p inputs does not contain the correct keys.
 *                            Strong throw guarantee.
 */
template<typename PropType>
auto unpack_basis_sets(const sde::type::input_map& inputs) {
    using element_type = double; // TODO: Get from PropType via TMP
    using basis_set    = integrals::type::basis_set<element_type>;
    using basis_vector = std::vector<basis_set>;
    using namespace property_types::ao_integrals;

    constexpr auto n_centers = n_centers_v<PropType>;
    auto spaces              = unpack_spaces<PropType>(inputs);

    basis_vector bs;
    for(const auto& space_i : spaces) bs.push_back(space_i.basis_set());
    return bs;
}

} // namespace integrals
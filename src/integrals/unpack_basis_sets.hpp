#pragma once
#include <property_types/ao_integrals/type_traits.hpp>
#include <property_types/ao_integrals/utilities/unpack_spaces.hpp>

namespace integrals {

/** @brief Unpacks the AO basis sets from the inptus to an integral module.
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
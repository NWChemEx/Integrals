#pragma once
#include "integrals/property_types.hpp"
#include "integrals/types.hpp"
#include <property_types/ao_integrals/type_traits.hpp>
#include <property_types/ao_integrals/utilities/unpack_spaces.hpp>
#include <sde/module_base.hpp>
#include <utilities/iter_tools.hpp>

namespace integrals {

template<typename T>
struct mode_shift : std::integral_constant<std::size_t, 0> {};

/// Specializes mode_shift so that multipoles skip the first mode
template<typename T, unsigned Order>
struct mode_shift<property_types::ao_integrals::EMultipole<T, Order>>
  : std::integral_constant<std::size_t, 1> {};

template<typename T>
static constexpr auto mode_shift_v = mode_shift<T>::value;

template<typename BaseType>
DECLARE_MODULE(Transformed);

template<typename BaseType>
TEMPLATED_MODULE_CTOR(Transformed, BaseType) {
    using my_pt = pt::transformed<BaseType>;

    satisfies_property_type<my_pt>();

    add_submodule<BaseType>("AO integral");
}

template<typename BaseType>
TEMPLATED_MODULE_RUN(Transformed, BaseType) {
    using my_pt        = pt::transformed<BaseType>;
    using element_type = double; // TODO: Get from TMP
    using ao_space_t   = type::ao_space_t<element_type>;
    using mo_space_t   = type::orbital_space_t<element_type>;
    using mo_ref_t     = std::reference_wrapper<const mo_space_t>;
    using mo_vec_t     = std::vector<mo_ref_t>;
    using namespace property_types::ao_integrals;

    // How many centers does this integral have?
    constexpr auto n_centers = n_centers_v<BaseType>;

    // Base property type for that many centers
    using n_center_t = NCenter<n_centers, element_type>;

    ////////////////////////////////////////////////////////////////////////////
    // Step 0: Get the integral we're going to transform
    ////////////////////////////////////////////////////////////////////////////

    // N.B. the inputs to the submod are a subset of the inputs to this mod, so
    // we just forward all the inputs...
    auto ao_results = submods.at("AO integral").value().run(inputs);
    auto [X]        = BaseType::unwrap_results(ao_results);

    ////////////////////////////////////////////////////////////////////////////
    // Step 1: Unpack the input and final spaces
    ////////////////////////////////////////////////////////////////////////////
    auto input_spaces = unpack_spaces<BaseType>(inputs);
    auto t            = my_pt::unwrap_inputs(inputs);
    mo_ref_t s0       = std::cref(std::get<0>(t));
    mo_ref_t s1       = std::cref(std::get<1>(t));
    mo_vec_t final_spaces{s0, s1}; // It has at least 2

    if constexpr(n_centers >= 3) {
        final_spaces.push_back(std::cref(std::get<2>(t)));
    }
    if constexpr(n_centers >= 4) {
        final_spaces.push_back(std::cref(std::get<3>(t)));
    }

    ////////////////////////////////////////////////////////////////////////////
    // Step 2: Determine which modes need transformed
    ////////////////////////////////////////////////////////////////////////////
    std::vector<bool> needs_transformed;
    for(const auto& [ii, fi] : utilities::Zip(input_spaces, final_spaces)) {
        // If the initial and final spaces aren't the same it needs
        // transformed
        needs_transformed.push_back(ii != fi.get());
    }

    ////////////////////////////////////////////////////////////////////////////
    // Step 3: Determine the order to do the transforms in
    ////////////////////////////////////////////////////////////////////////////

    /* Our strategy here is to do the transformations in the order that
     * reduces the size the quickest (or expands the size the slowest).
     *
     * Let  d = final_size - initial_size, then a positive d means we're
     * enlarging, negative d means we're reducing.
     *
     * A multimap with d's as keys and modes as values will tell us the order to
     * do the transformations in (this relies on the fact that multimap stores
     * elements sorted by key).
     */

    std::multimap<long int, std::size_t> order;
    for(auto i = 0; i < n_centers; ++i) {
        if(!needs_transformed[i]) continue;
        const auto final_size = final_spaces[i].get().size();
        const auto init_size  = input_spaces[i].size();
        long int d            = 0;
        if(final_size < init_size) {
            d = init_size - final_size;
            d *= -1;
        } else {
            d = final_size - init_size;
        }
        order.emplace(std::make_pair(d, i));
    }
    // Step 3: Transform the tensor

    constexpr auto shift = mode_shift_v<BaseType>;
    for(auto [d, mode] : order) {
        X = final_spaces[mode].get().transform(X, std::vector{mode + shift});
    }

    auto rv = results();
    return my_pt::wrap_results(rv, X);
}

/** @brief Wraps the process of registering a transformed integral.
 *
 *  The Transformed module is templated on the property type being transformed.
 *  This means that a new module needs to be created for each property type.
 *  This function will instantiate that module and add it to the given module
 *  manager under the key `Transformed X` where `X` is the key of the module
 *  being called to generate the AO integral. Strictly speaking the same
 *  `Transformed<X>` can be used for all modules satisfying `X`, but it's just
 *  easier at the moment to make one for each module to avoid having to switch
 *  submodules.
 *
 *  @param[in,out] mm The ModuleManager instance we are adding the new module
 to.
 *                 The submodule to use for building the AO integrals must also
 *                 live in @p mm. After this call @p mm will contain the newly
                   created module.
 *  @param[in] ao_key The key of the module to use for building the AO
               integrals. Must live in @p mm.
 */
template<typename PropType>
void register_transformed_integral(sde::ModuleManager& mm,
                                   const std::string& ao_key) {
    const auto new_key = std::string("Transformed ") + ao_key;
    mm.add_module<Transformed<PropType>>(new_key);
    mm.change_submod(new_key, "AO integral", ao_key);
}

extern template class Transformed<pt::correlation_factor_4c<double>>;
extern template class Transformed<pt::doi<double>>;
extern template class Transformed<pt::eri2c<double>>;

} // namespace integrals
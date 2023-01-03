/*
 * Copyright 2023 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * Copyright 2022 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "integrals/libint_integral.hpp"
#include "nuclear_integral.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include <property_types/ao_integrals/nuclear.hpp>

namespace integrals {

template<typename element_type>
using nuclear_type = property_types::ao_integrals::Nuclear<element_type>;
template<typename element_type>
using libint_type = property_types::LibIntIntegral<element_type>;
template<typename element_type>
using tensor = typename type::tensor<element_type>;
template<typename element_type>
using value_type = typename tensor<element_type>::value_type;

template<typename element_type>
NuclearInt<element_type>::NuclearInt() : sde::ModuleBase(this) {
    description("Computes nuclear integrals with Libint");
    satisfies_property_type<nuclear_type<element_type>>();
    satisfies_property_type<libint_type<element_type>>();
}

template<typename element_type>
sde::type::result_map NuclearInt<element_type>::run_(
  sde::type::input_map inputs, sde::type::submodule_map submods) const {
    auto [mol, bra_space, ket_space] =
      nuclear_type<element_type>::unwrap_inputs(inputs);

    auto& bra         = bra_space.basis_set();
    auto& ket         = ket_space.basis_set();
    std::size_t deriv = 0;

    auto [thresh, tile_size, cs_thresh, atom_ranges] =
      libint_type<element_type>::unwrap_inputs(inputs);
    auto& world = TA::get_default_world();

    std::vector<std::pair<double, std::array<double, 3>>> qs;
    for(const auto& ai : mol)
        qs.emplace_back(static_cast<const double&>(ai.Z()), ai.coords());

    auto fill = nwx_TA::FillNDFunctor<value_type<element_type>,
                                      libint2::Operator::nuclear, 2>();
    fill.initialize(nwx_libint::make_basis_sets({bra, ket}), deriv, thresh,
                    cs_thresh);
    fill.factory.qs = qs;

    auto trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges);

    auto V = TiledArray::make_array<tensor<element_type>>(world, trange, fill);

    auto rv = results();
    return nuclear_type<element_type>::wrap_results(rv, V);
}

template class NuclearInt<double>;

} // namespace integrals

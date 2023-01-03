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

#include "../unpack_basis_sets.hpp"
#include "fill_ND_functor.hpp"
#include "libint.hpp"
#include "nwx_TA_utils.hpp"
#include "nwx_libint.hpp"
#include "traits.hpp"
#include <chemist/ta_helpers/ta_hashers.hpp>
#include <property_types/ao_integrals/type_traits.hpp>
#include <stdexcept>

namespace integrals {

template<typename PropType>
TEMPLATED_MODULE_CTOR(Libint, PropType) {
    using element_type = double; // TODO: Get from PropType
    using type::pair_vector;
    using type::size_vector;

    description("Computes an in-core integral with libint");
    satisfies_property_type<PropType>();

    add_input<element_type>("Threshold").set_default(1.0E-16);
    add_input<size_vector>("Tile Size").set_default(size_vector{180});
    add_input<pair_vector>("Atom Tile Groups").set_default(pair_vector{});
}

template<typename PropType>
TEMPLATED_MODULE_RUN(Libint, PropType) {
    using element_type = double; // TODO: Get from PropType
    using tensor_type  = type::tensor<element_type>;
    using value_type   = typename tensor_type::value_type;
    using type::pair_vector;
    using type::size_vector;

    auto& world      = TA::get_default_world(); // TODO: Get from runtime
    auto thresh      = inputs.at("Threshold").value<element_type>();
    auto tile_size   = inputs.at("Tile Size").value<size_vector>();
    auto atom_ranges = inputs.at("Atom Tile Groups").value<pair_vector>();

    constexpr auto n_centers =
      property_types::ao_integrals::n_centers_v<PropType>;
    constexpr auto op = op_v<PropType>;

    auto bs = unpack_basis_sets<PropType>(inputs);
    if constexpr(property_types::ao_integrals::is_doi_v<PropType>) {
        decltype(bs){bs[0], bs[0], bs[1], bs[1]}.swap(bs);
    }
    auto trange = nwx_TA::select_tiling(bs, tile_size, atom_ranges);
    auto fill   = nwx_TA::FillNDFunctor<value_type, op, n_centers>();
    const std::size_t deriv      = 0;   // TODO: Template on derivative order
    const element_type cs_thresh = 0.0; // Just to satisfy initialize
    fill.initialize(nwx_libint::make_basis_sets(bs), deriv, thresh, cs_thresh);

    // Take care of any special parameters the fill function needs
    // (should probably we done inside FillNDFunctor to encapsulate the setup)
    if constexpr(property_types::ao_integrals::is_nuclear_v<PropType>) {
        using mol_type  = const chemist::Molecule&;
        const auto& mol = inputs.at("Molecule").value<mol_type>();
        std::vector<std::pair<double, std::array<double, 3>>> qs;
        for(const auto& ai : mol)
            qs.emplace_back(static_cast<const double&>(ai.Z()), ai.coords());
        fill.factory.qs = qs;
    } else if constexpr(property_types::ao_integrals::is_stg_v<PropType> ||
                        property_types::ao_integrals::is_yukawa_v<PropType>) {
        auto gamma = inputs.at("STG Exponent").value<element_type>();
        fill.factory.stg_exponent = gamma;
    }

    auto I  = TiledArray::make_array<tensor_type>(world, trange, fill);
    auto rv = results();
    return PropType::wrap_results(rv, I);
}

template class Libint<pt::doi<double>>;
template class Libint<pt::eri2c<double>>;
template class Libint<pt::eri3c<double>>;
template class Libint<pt::eri4c<double>>;
template class Libint<pt::kinetic<double>>;
template class Libint<pt::nuclear<double>>;
template class Libint<pt::overlap<double>>;
template class Libint<pt::stg2c<double>>;
template class Libint<pt::stg3c<double>>;
template class Libint<pt::stg4c<double>>;
template class Libint<pt::yukawa2c<double>>;
template class Libint<pt::yukawa3c<double>>;
template class Libint<pt::yukawa4c<double>>;

} // namespace integrals

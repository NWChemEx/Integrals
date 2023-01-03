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

#include "integrals/eri_integrals.hpp"
#include "integrals/libint_integral.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include <property_types/ao_integrals/electron_repulsion.hpp>
#include <property_types/cauchy_schwarz_approximation.hpp>

namespace integrals {

template<typename element_type>
using eri3c_type = property_types::ao_integrals::ERI3C<element_type>;
template<typename element_type>
using libint_type = property_types::LibIntIntegral<element_type>;
template<typename element_type>
using tensor = typename type::tensor<element_type>;
template<typename element_type>
using value_type = typename tensor<element_type>::value_type;
template<typename element_type>
using cs_approx_type = property_types::CauchySchwarzApprox<element_type>;

template<typename element_type>
ERI3CInt<element_type>::ERI3CInt() : sde::ModuleBase(this) {
    description("Computes 3-center electron repulsion integrals with Libint");
    satisfies_property_type<eri3c_type<element_type>>();
    satisfies_property_type<libint_type<element_type>>();

    add_submodule<cs_approx_type<element_type>>("Cauchy-Schwarz")
      .set_description(
        "Computes the Cauchy-Schwarz Matrix for a pair of basis sets");
}

template<typename element_type>
sde::type::result_map ERI3CInt<element_type>::run_(
  sde::type::input_map inputs, sde::type::submodule_map submods) const {
    auto [bra_space, ket1_space, ket2_space] =
      eri3c_type<element_type>::unwrap_inputs(inputs);

    auto& bra         = bra_space.basis_set();
    auto& ket1        = ket1_space.basis_set();
    auto& ket2        = ket2_space.basis_set();
    std::size_t deriv = 0;

    auto [thresh, tile_size, cs_thresh, atom_ranges] =
      libint_type<element_type>::unwrap_inputs(inputs);
    auto& world = TA::get_default_world();

    auto fill = nwx_TA::FillNDFunctor<value_type<element_type>,
                                      libint2::Operator::coulomb, 3>();
    fill.initialize(nwx_libint::make_basis_sets({bra, ket1, ket2}), deriv,
                    thresh, cs_thresh);

    if(cs_thresh > 0.0) {
        auto [cs_mat] =
          submods.at("Cauchy-Schwarz")
            .run_as<cs_approx_type<element_type>>(ket1, ket2, deriv);
        fill.screen.cs_mat2 = cs_mat;
    }

    auto trange =
      nwx_TA::select_tiling({bra, ket1, ket2}, tile_size, atom_ranges);

    auto I = TiledArray::make_array<tensor<element_type>>(world, trange, fill);

    auto rv = results();
    return eri3c_type<element_type>::wrap_results(rv, I);
}

template class ERI3CInt<double>;

} // namespace integrals

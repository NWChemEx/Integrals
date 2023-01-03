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

#include "integrals/emultipole_integrals.hpp"
#include "integrals/libint_integral.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include <chemist/ta_helpers/einsum/einsum.hpp>
#include <property_types/ao_integrals/emultipole.hpp>
#include <property_types/ao_integrals/overlap.hpp>

namespace integrals {

template<typename element_type>
using overlap_type = property_types::ao_integrals::Overlap<element_type>;
template<typename element_type>
using eDipole_type = property_types::ao_integrals::EDipole<element_type>;
template<typename element_type>
using eQuadrupole_type =
  property_types::ao_integrals::EQuadrupole<element_type>;
template<typename element_type>
using eOctopole_type = property_types::ao_integrals::EOctopole<element_type>;
template<typename element_type>
using libint_type = property_types::LibIntIntegral<element_type>;
template<typename element_type>
using tensor = typename type::tensor<element_type>;
template<typename element_type>
using value_type = typename tensor<element_type>::value_type;

template<typename element_type>
EOctopoleInt<element_type>::EOctopoleInt() : sde::ModuleBase(this) {
    description("Computes dipole integrals with Libint");
    satisfies_property_type<overlap_type<element_type>>();
    satisfies_property_type<eDipole_type<element_type>>();
    satisfies_property_type<eQuadrupole_type<element_type>>();
    satisfies_property_type<eOctopole_type<element_type>>();
    satisfies_property_type<libint_type<element_type>>();
}

template<typename element_type>
sde::type::result_map EOctopoleInt<element_type>::run_(
  sde::type::input_map inputs, sde::type::submodule_map submods) const {}

template class EOctopoleInt<double>;

} // namespace integrals

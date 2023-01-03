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

#pragma once
#include <simde/simde.hpp>

namespace integrals {

template<std::size_t N, typename OpType>
DECLARE_MODULE(StandardTransform);

extern template class StandardTransform<2, simde::type::el_scf_k>;
extern template class StandardTransform<2, simde::type::fock>;
extern template class StandardTransform<2, simde::type::el_kinetic>;
extern template class StandardTransform<2, simde::type::el_nuc_coulomb>;
extern template class StandardTransform<2, simde::type::fock>;
extern template class StandardTransform<3, simde::type::el_el_coulomb>;
extern template class StandardTransform<4, simde::type::el_el_coulomb>;
extern template class StandardTransform<4, simde::type::el_el_f12_commutator>;
extern template class StandardTransform<4, simde::type::el_el_stg>;
extern template class StandardTransform<4, simde::type::el_el_yukawa>;
} // namespace integrals

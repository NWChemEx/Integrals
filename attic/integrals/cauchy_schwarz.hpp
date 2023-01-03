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
#include <libint2.hpp>
#include <sde/module_base.hpp>

namespace integrals {

template<typename element_type = double,
         libint2::Operator op  = libint2::Operator::coulomb>
DECLARE_MODULE(CauchySchwarz);

extern template class CauchySchwarz<double>;
extern template class CauchySchwarz<double, libint2::Operator::stg>;
extern template class CauchySchwarz<double, libint2::Operator::yukawa>;

using CS_ERI    = CauchySchwarz<double>;
using CS_STG    = CauchySchwarz<double, libint2::Operator::stg>;
using CS_Yukawa = CauchySchwarz<double, libint2::Operator::yukawa>;

} // namespace integrals
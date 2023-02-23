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

#include "ao_integrals/ao_integrals.hpp"
#include "libint/libint.hpp"
#include "transforms/transforms.hpp"
#include <integrals/integrals_mm.hpp>

using namespace simde::type;

namespace integrals {

void load_modules(pluginplay::ModuleManager& mm) {
    ao_integrals::load_ao_integrals(mm);
    libint::load_libint_modules(mm);
    transforms::load_transformed_libint_integrals(mm);

    ao_integrals::ao_integrals_set_defaults(mm);
}

} // namespace integrals

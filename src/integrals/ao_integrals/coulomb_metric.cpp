/*
 * Copyright 2025 NWChemEx-Project
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

#include "ao_integrals.hpp"

namespace integrals::ao_integrals {

using pt = simde::ERI2;

namespace {

auto desc = R"(
Inverse Coulomb Metric
---------------------
)";

}

MODULE_CTOR(CoulombMetric) {
    description(desc);
    satisfies_property_type<pt>();
    add_submodule<pt>("Two-center ERI");
}

MODULE_RUN(CoulombMetric) {
    const auto& [braket] = pt::unwrap_inputs(inputs);
    auto& eri2_mod       = submods.at("Two-center ERI");

    const auto& M = eri2_mod.run_as<pt>(braket);

    // Cholesky Decomp

    auto rv = results();
    return pt::wrap_results(rv, M);
}

} // namespace integrals::ao_integrals
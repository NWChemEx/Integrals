/*
 * Copyright 2026 NWChemEx-Project
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
#include <integrals/integrals.hpp>

namespace integrals::ao_integrals {
namespace {

const auto desc = R"(
UQ Initializer
--------------

Produces a UQFactory configured by this module's "UQ Type" input. The
resulting factory can be used directly in a module's hot loops, avoiding
module-dispatch overhead.

The modules which create uncertain values from double values, e.g., the
one-electron integrals, share this module, so that the kind of uncertainty is
set in one place. The default, "none", makes double values.
)";

} // namespace

using pt = simde::UQInitializer;

MODULE_CTOR(UQInitializer) {
    satisfies_property_type<pt>();
    description(desc);
    add_input<std::string>("UQ Type").set_default("none").set_description(
      "The kind of uncertainty: none, uncertain, interval, "
      "affine, or thresholded affine");
}

MODULE_RUN(UQInitializer) {
    auto uq_type = inputs.at("UQ Type").value<std::string>();
    auto kind    = simde::uq_kind_from_string(uq_type);
    simde::UQFactory factory(kind);
    auto rv = results();
    return pt::wrap_results(rv, factory);
}

} // namespace integrals::ao_integrals

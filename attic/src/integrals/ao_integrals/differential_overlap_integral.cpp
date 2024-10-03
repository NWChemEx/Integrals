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

#include "ao_integrals.hpp"
#include "detail_/aos2shells.hpp"
#include "detail_/get_coeff.hpp"
#include "detail_/make_shape.hpp"
#include "detail_/select_allocator.hpp"
#include "detail_/shells2ord.hpp"
#include "detail_/unpack_bases.hpp"
#include <simde/integral_factory.hpp>
#include <simde/tensor_representation/ao_tensor_representation.hpp>

namespace integrals::ao_integrals {

/// Grab the various detail_ functions
using namespace detail_;

MODULE_CTOR(AOIntegralDOI) {
    description("Computes DOI integrals");
    using op_t    = simde::type::el_el_delta;
    using my_pt   = simde::AOTensorRepresentation<2, op_t>;
    using doi4_pt = simde::AOTensorRepresentation<4, op_t>;

    satisfies_property_type<my_pt>();

    add_submodule<doi4_pt>("DOI4").set_description("Computes DOI4 integrals");
}

MODULE_RUN(AOIntegralDOI) {
    using op_t    = simde::type::el_el_delta;
    using my_pt   = simde::AOTensorRepresentation<2, op_t>;
    using doi4_pt = simde::AOTensorRepresentation<4, op_t>;

    /// Get inputs
    const auto& [bra, op, ket] = my_pt::unwrap_inputs(inputs);
    auto& doi4_mod             = submods.at("DOI4");

    /// Run DOI4
    auto I = doi4_mod.run_as<doi4_pt>(bra, bra, op, ket, ket);

    /// Finish
    auto rv = results();
    return my_pt::wrap_results(rv, I);
}

} // namespace integrals::ao_integrals

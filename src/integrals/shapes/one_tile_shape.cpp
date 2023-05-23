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

#include "shapes.hpp"
#include <simde/simde.hpp>

namespace integrals::shapes {

using integral_shape_t = simde::IntegralShape;
using shape_t          = typename simde::type::tensor::shape_type;
using extents_t        = typename shape_t::extents_type;

MODULE_CTOR(OneTileShape) {
    satisfies_property_type<integral_shape_t>();
    description("Construct a tensor shape that has one big tile");
}

MODULE_RUN(OneTileShape) {
    auto [bases] = integral_shape_t::unwrap_inputs(inputs);

    extents_t extents;
    for(auto& set : bases) { extents.push_back(set.n_aos()); }
    shape_t shape{extents};

    auto rv = results();
    return integral_shape_t::wrap_results(rv, shape);
}

} // namespace integrals::shapes

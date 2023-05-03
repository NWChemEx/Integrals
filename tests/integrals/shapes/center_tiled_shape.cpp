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

#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <integrals/property_types/integral_shape.hpp>
#include <mokup/mokup.hpp>

using pt      = integrals::IntegralShape;
using shape_t = typename simde::type::tensor::shape_type;
using input_t = std::vector<simde::type::ao_basis_set>;

using namespace mokup;

TEST_CASE("CenterTiledShape") {
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    auto name         = molecule::h2o;
    auto sto3g_space  = get_bases(name, basis_set::sto3g);
    auto ccpvdz_space = get_bases(name, basis_set::ccpvdz);

    input_t inputs(2);
    inputs[0] = sto3g_space.basis_set();
    inputs[1] = ccpvdz_space.basis_set();

    auto shape = mm.at("CenterTiledShape").run_as<pt>(inputs);
    shape_t corr{{{0, 5, 6, 7}, {0, 14, 19, 24}}};
    REQUIRE(shape == corr);
}

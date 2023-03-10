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

#include "direct_allclose.hpp"
#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>
#include <simde/tensor_representation/tensor_representation.hpp>
#include <tensorwrapper/tensor/allclose.hpp>

using namespace mokup;

TEST_CASE("Yukawa2C") {
    using op_type       = simde::type::el_el_yukawa;
    using integral_type = simde::AOTensorRepresentation<2, op_type>;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    auto name = molecule::h2o;
    auto bs   = basis_set::sto3g;
    auto aos  = get_bases(name, bs);
    std::vector bases{bs, bs};
    auto corr = get_ao_data(name, bases, property::yukawa);

    chemist::Electron e;
    op_type gr(chemist::operators::STG(1.0, 1.0), e, e);

    SECTION("Explicit") {
        auto [X] = mm.at("Yukawa2").run_as<integral_type>(aos, gr, aos);
        REQUIRE(tensorwrapper::tensor::allclose(X, corr));
    }

    SECTION("Direct") {
        auto [X] = mm.at("Direct Yukawa2").run_as<integral_type>(aos, gr, aos);
        REQUIRE(direct_allclose(X, corr));
    }
}
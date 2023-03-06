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

TEST_CASE("Overlap") {
    using op_type       = simde::type::el_identity;
    using integral_type = simde::AOTensorRepresentation<2, op_type>;

    // TODO: Better names which describe what's being tested
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto aos        = get_bases(name, bs);
    std::vector bases{bs, bs};
    auto corr_S = get_ao_data(name, bases, property::overlap);

    op_type I;

    SECTION("Explicit") {
        auto [S] = mm.at("Overlap").run_as<integral_type>(aos, I, aos);
        REQUIRE(tensorwrapper::tensor::allclose(S, corr_S));
    }

    SECTION("Direct") {
        auto [S] = mm.at("Direct Overlap").run_as<integral_type>(aos, I, aos);
        REQUIRE(direct_allclose(S, corr_S));
    }
}

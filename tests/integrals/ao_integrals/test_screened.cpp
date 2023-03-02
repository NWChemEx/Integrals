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
#include <simde/simde.hpp>
#include <tensorwrapper/tensor/allclose.hpp>

using namespace mokup;

TEST_CASE("Overlap CS") {
    using integral_type = simde::EOverlap;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    mm.change_input("Overlap CS", "Screening Threshold", 0.5);
    mm.change_input("Direct Overlap CS", "Screening Threshold", 0.5);

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto aos        = get_bases(name, bs);
    std::vector bases{bs, bs};
    auto corr_S = get_ao_data(name, bases, property::screened_overlap);

    simde::type::el_identity I;

    SECTION("Explicit") {
        auto [S] = mm.at("Overlap CS").run_as<integral_type>(aos, I, aos);
        REQUIRE(tensorwrapper::tensor::allclose(S, corr_S));
    }

    SECTION("Direct") {
        auto [S] =
          mm.at("Direct Overlap CS").run_as<integral_type>(aos, I, aos);
        REQUIRE(direct_allclose(S, corr_S));
    }
}

TEST_CASE("ERI4C CS") {
    using integral_type = simde::ERI4;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    mm.change_input("ERI4 CS", "Screening Threshold", 0.005);
    mm.change_input("Direct ERI4 CS", "Screening Threshold", 0.005);

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto aos        = get_bases(name, bs);
    std::vector bases{bs, bs, bs, bs};
    auto corr_X = get_ao_data(name, bases, property::screened_eris);

    simde::type::el_el_coulomb r12;

    SECTION("Explicit") {
        auto [X] =
          mm.at("ERI4 CS").run_as<integral_type>(aos, aos, r12, aos, aos);
        REQUIRE(tensorwrapper::tensor::allclose(X, corr_X));
    }

    SECTION("Direct") {
        auto [X] = mm.at("Direct ERI4 CS")
                     .run_as<integral_type>(aos, aos, r12, aos, aos);
        REQUIRE(direct_allclose(X, corr_X));
    }
}

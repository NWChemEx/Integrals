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
#include <mokup/mokup.hpp>
#include <simde/tensor_representation/tensor_representation.hpp>
#include <tensorwrapper/tensor/allclose.hpp>

using namespace integrals;
using namespace mokup;

TEST_CASE("Quadrupole") {
    using i_op   = simde::type::el_identity;
    using d_op   = simde::type::el_dipole;
    using q_op   = simde::type::el_quadrupole;
    using s_type = simde::AOTensorRepresentation<2, i_op>;
    using d_type = simde::AOTensorRepresentation<2, d_op>;
    using q_type = simde::AOTensorRepresentation<2, q_op>;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto aos        = get_bases(name, bs);
    std::vector bases{bs, bs};
    d_op r;
    q_op r2;

    /// TODO: this needs to actually test something.
    // SECTION("overlap matrix") {
    //     mm.change_input("EDipole", "Origin", origin);
    //     auto [S] = mm.at("EDipole").run_as<s_type>(aos, aos);
    //     REQUIRE(tensorwrapper::ta_helpers::allclose(S, X));
    // }

    SECTION("dipole matrix") {
        auto D  = mm.at("EQuadrupole").run_as<d_type>(aos, r, aos);
        auto corr = get_ao_data(name, bases, property::dipole);
        REQUIRE(tensorwrapper::tensor::allclose(D, corr));
    }

    SECTION("Quadrupole") {
        auto Q  = mm.at("EQuadrupole").run_as<q_type>(aos, r2, aos);
        auto corr = get_ao_data(name, bases, property::quadrupole);
        REQUIRE(tensorwrapper::tensor::allclose(Q, corr));
    }
}

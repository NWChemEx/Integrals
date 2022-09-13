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
#include <tensorwrapper/tensor/allclose.hpp>

using namespace mokup;

TEST_CASE("ERI2C") {
    using op            = simde::type::el_el_coulomb;
    using integral_type = simde::AOTensorRepresentation<2, op>;
    using size_vector   = std::vector<std::size_t>;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    auto name = molecule::h2o;
    auto bs   = basis_set::sto3g;
    auto mol  = get_molecule(name);
    auto aos  = get_bases(name, bs);
    std::vector bases{bs, bs};
    auto corr = get_ao_data(name, bases, property::eris);

    // mm.at("ERI2").change_input("Tile size", size_vector{1});
    simde::type::el_el_coulomb r12;
    auto [X] = mm.at("ERI2").run_as<integral_type>(aos, r12, aos);
    REQUIRE(tensorwrapper::tensor::allclose(X, corr));
}

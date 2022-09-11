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

TEST_CASE("Kinetic") {
    using op_type       = simde::type::el_kinetic;
    using integral_type = simde::AOTensorRepresentation<2, op_type>;
    using size_vector   = std::vector<std::size_t>;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    auto name = molecule::h2o;
    auto bs   = basis_set::sto3g;
    auto aos  = get_bases(name, bs);
    std::vector bases{bs, bs};
    auto corr = get_ao_data(name, bases, property::kinetic);

    // mm.at("Kinetic").change_input("Tile size", size_vector{1, 2});
    op_type t;
    auto [T] = mm.at("Kinetic").run_as<integral_type>(aos, t, aos);
    REQUIRE(tensorwrapper::tensor::allclose(T, corr));
}

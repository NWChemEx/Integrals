// /*
//  * Copyright 2022 NWChemEx-Project
//  *
//  * Licensed under the Apache License, Version 2.0 (the "License");
//  * you may not use this file except in compliance with the License.
//  * You may obtain a copy of the License at
//  *
//  * http://www.apache.org/licenses/LICENSE-2.0
//  *
//  * Unless required by applicable law or agreed to in writing, software
//  * distributed under the License is distributed on an "AS IS" BASIS,
//  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  * See the License for the specific language governing permissions and
//  * limitations under the License.
//  */

// #include "integrals/integrals.hpp"
// #include <catch2/catch.hpp>
// #include <mokup/mokup.hpp>
// #include <simde/tensor_representation/tensor_representation.hpp>
// #include <tensorwrapper/tensor/allclose.hpp>

// using namespace mokup;

// TEST_CASE("STG 2 Center Correlation Factor") {
//     using op_type       = simde::type::el_el_stg;
//     using integral_type = simde::AOTensorRepresentation<2, op_type>;
//     const auto key      = "STG2";

//     pluginplay::ModuleManager mm;
//     integrals::load_modules(mm);

//     const auto name = molecule::h2;
//     const auto prop = property::stg_correlation_factor;

//     chemist::Electron e;
//     op_type stg;

//     for(const auto& bs : {basis_set::sto3g, basis_set::ccpvdz}) {
//         std::vector<mokup::basis_set> bs_key(2, bs);
//         SECTION(as_string(name, bs)) {
//             auto aos    = get_bases(name, bs);
//             auto X_corr = get_ao_data(name, bs_key, prop);
//             auto [X]    = mm.at(key).run_as<integral_type>(aos, stg, aos);
//             REQUIRE(tensorwrapper::tensor::allclose(X, X_corr));
//         }
//     }
// }

// TEST_CASE("STG 4 Center Correlation Factor") {
//     using op_type       = simde::type::el_el_stg;
//     using integral_type = simde::AOTensorRepresentation<4, op_type>;
//     const auto key      = "STG4";

//     pluginplay::ModuleManager mm;
//     integrals::load_modules(mm);

//     const auto name = molecule::h2;
//     const auto prop = property::stg_correlation_factor;
//     chemist::Electron e;
//     op_type stg;

//     for(const auto& bs : {basis_set::sto3g, basis_set::ccpvdz}) {
//         std::vector<basis_set> bs_key(4, bs);
//         SECTION(as_string(name, bs)) {
//             auto aos    = get_bases(name, bs);
//             auto X_corr = get_ao_data(name, bs_key, prop);
//             auto [X] =
//               mm.at(key).run_as<integral_type>(aos, aos, stg, aos, aos);
//             REQUIRE(tensorwrapper::tensor::allclose(X, X_corr));
//         }
//     }
// }

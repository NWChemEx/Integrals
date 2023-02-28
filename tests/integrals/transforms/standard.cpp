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
// #include <simde/simde.hpp>
// #include <tensorwrapper/tensor/allclose.hpp>

// using namespace mokup;

// TEST_CASE("Transformed") {
//     using op_type  = simde::type::el_el_coulomb;
//     using submod_t = simde::AOTensorRepresentation<4, op_type>;
//     using pt       = simde::TransformedTensorRepresentation<4, op_type>;

//     using ao_map = simde::space_map_t<simde::type::ao_space>;
//     using mo_map = simde::space_map_t<simde::type::derived_space>;

//     pluginplay::ModuleManager mm;
//     integrals::load_modules(mm);

//     auto name = molecule::h2;
//     auto bs   = basis_set::sto3g;
//     auto aos  = get_bases(name, bs);
//     auto mos  = get_space(property::occupied, name, bs);
//     std::vector bases{bs, bs, bs, bs};
//     auto G = get_ao_data(name, bases, property::eris);

//     auto& mod     = mm.at("Transformed ERI4");
//     const auto& C = mos.C();
//     op_type r12;

//     auto submod = pluginplay::make_lambda<submod_t>(
//       [&](auto&& b1, auto&& b2, auto&& op, auto&& k1, auto&& k2) {
//           REQUIRE(b1 == aos);
//           REQUIRE(b2 == aos);
//           REQUIRE(op == r12);
//           REQUIRE(k1 == aos);
//           REQUIRE(k2 == aos);
//           return G;
//       });
//     mod.change_submod("Integral Kernel", submod);

//     SECTION("Transform mode 0") {
//         ao_map ao_spaces;
//         mo_map mo_spaces;
//         mo_spaces.emplace(0, mos);
//         ao_spaces.emplace(1, aos);
//         ao_spaces.emplace(2, aos);
//         ao_spaces.emplace(3, aos);
//         auto [rv] = mod.run_as<pt>(ao_spaces, mo_spaces, r12);
//         simde::type::tensor corr;
//         corr("i,b,c,d") = C("a,i") * G("a,b,c,d");
//         REQUIRE(tensorwrapper::tensor::allclose(rv, corr));
//     }

//     SECTION("Transform mode 1") {
//         ao_map ao_spaces;
//         mo_map mo_spaces;
//         mo_spaces.emplace(1, mos);
//         ao_spaces.emplace(0, aos);
//         ao_spaces.emplace(2, aos);
//         ao_spaces.emplace(3, aos);
//         auto [rv] = mod.run_as<pt>(ao_spaces, mo_spaces, r12);
//         simde::type::tensor corr;
//         corr("a,i,c,d") = C("b,i") * G("a,b,c,d");
//         REQUIRE(tensorwrapper::tensor::allclose(rv, corr));
//     }

//     SECTION("Transform mode 2") {
//         ao_map ao_spaces;
//         mo_map mo_spaces;
//         mo_spaces.emplace(2, mos);
//         ao_spaces.emplace(0, aos);
//         ao_spaces.emplace(1, aos);
//         ao_spaces.emplace(3, aos);
//         auto [rv] = mod.run_as<pt>(ao_spaces, mo_spaces, r12);
//         simde::type::tensor corr;
//         corr("a,b,i,d") = C("c,i") * G("a,b,c,d");
//         REQUIRE(tensorwrapper::tensor::allclose(rv, corr));
//     }

//     SECTION("Transform mode 3") {
//         ao_map ao_spaces;
//         mo_map mo_spaces;
//         mo_spaces.emplace(3, mos);
//         ao_spaces.emplace(0, aos);
//         ao_spaces.emplace(1, aos);
//         ao_spaces.emplace(2, aos);
//         auto [rv] = mod.run_as<pt>(ao_spaces, mo_spaces, r12);
//         simde::type::tensor corr;
//         corr("a,b,c,i") = C("d,i") * G("a,b,c,d");
//         REQUIRE(tensorwrapper::tensor::allclose(rv, corr));
//     }

//     SECTION("Transform modes 0 and 1") {
//         ao_map ao_spaces;
//         mo_map mo_spaces;
//         mo_spaces.emplace(0, mos);
//         mo_spaces.emplace(1, mos);
//         ao_spaces.emplace(2, aos);
//         ao_spaces.emplace(3, aos);
//         auto [rv] = mod.run_as<pt>(ao_spaces, mo_spaces, r12);
//         simde::type::tensor temp, corr;
//         temp("i,b,c,d") = C("a,i") * G("a,b,c,d");
//         corr("i,j,c,d") = C("b,j") * temp("i,b,c,d");
//         REQUIRE(tensorwrapper::tensor::allclose(rv, corr));
//     }

//     SECTION("Transform modes 0 and 2") {
//         ao_map ao_spaces;
//         mo_map mo_spaces;
//         mo_spaces.emplace(0, mos);
//         mo_spaces.emplace(2, mos);
//         ao_spaces.emplace(1, aos);
//         ao_spaces.emplace(3, aos);
//         auto [rv] = mod.run_as<pt>(ao_spaces, mo_spaces, r12);
//         simde::type::tensor temp, corr;
//         temp("i,b,c,d") = C("a,i") * G("a,b,c,d");
//         corr("i,b,j,d") = C("c,j") * temp("i,b,c,d");
//         REQUIRE(tensorwrapper::tensor::allclose(rv, corr));
//     }

//     SECTION("Transform modes 0 and 3") {
//         ao_map ao_spaces;
//         mo_map mo_spaces;
//         mo_spaces.emplace(0, mos);
//         mo_spaces.emplace(3, mos);
//         ao_spaces.emplace(1, aos);
//         ao_spaces.emplace(2, aos);
//         auto [rv] = mod.run_as<pt>(ao_spaces, mo_spaces, r12);
//         simde::type::tensor temp, corr;
//         temp("i,b,c,d") = C("a,i") * G("a,b,c,d");
//         corr("i,b,c,j") = C("d,j") * temp("i,b,c,d");
//         REQUIRE(tensorwrapper::tensor::allclose(rv, corr));
//     }

//     SECTION("Transform modes 1 and 2") {
//         ao_map ao_spaces;
//         mo_map mo_spaces;
//         mo_spaces.emplace(1, mos);
//         mo_spaces.emplace(2, mos);
//         ao_spaces.emplace(0, aos);
//         ao_spaces.emplace(3, aos);
//         auto [rv] = mod.run_as<pt>(ao_spaces, mo_spaces, r12);
//         simde::type::tensor temp, corr;
//         temp("a,i,c,d") = C("b,i") * G("a,b,c,d");
//         corr("a,i,j,d") = C("c,j") * temp("a,i,c,d");
//         REQUIRE(tensorwrapper::tensor::allclose(rv, corr));
//     }

//     SECTION("Transform modes 1 and 3") {
//         ao_map ao_spaces;
//         mo_map mo_spaces;
//         mo_spaces.emplace(1, mos);
//         mo_spaces.emplace(3, mos);
//         ao_spaces.emplace(0, aos);
//         ao_spaces.emplace(2, aos);
//         auto [rv] = mod.run_as<pt>(ao_spaces, mo_spaces, r12);
//         simde::type::tensor temp, corr;
//         temp("a,i,c,d") = C("b,i") * G("a,b,c,d");
//         corr("a,i,c,j") = C("d,j") * temp("a,i,c,d");
//         REQUIRE(tensorwrapper::tensor::allclose(rv, corr));
//     }

//     SECTION("Transform modes 2 and 3") {
//         ao_map ao_spaces;
//         mo_map mo_spaces;
//         mo_spaces.emplace(2, mos);
//         mo_spaces.emplace(3, mos);
//         ao_spaces.emplace(0, aos);
//         ao_spaces.emplace(1, aos);
//         auto [rv] = mod.run_as<pt>(ao_spaces, mo_spaces, r12);
//         simde::type::tensor temp, corr;
//         temp("a,b,i,d") = C("c,i") * G("a,b,c,d");
//         corr("a,b,i,j") = C("d,j") * temp("a,b,i,d");
//         REQUIRE(tensorwrapper::tensor::allclose(rv, corr));
//     }

//     SECTION("Transform modes 0, 1, and 2") {
//         ao_map ao_spaces;
//         mo_map mo_spaces;
//         mo_spaces.emplace(0, mos);
//         mo_spaces.emplace(1, mos);
//         mo_spaces.emplace(2, mos);
//         ao_spaces.emplace(3, aos);
//         auto [rv] = mod.run_as<pt>(ao_spaces, mo_spaces, r12);
//         simde::type::tensor temp, corr;
//         corr("i,b,c,d") = C("a,i") * G("a,b,c,d");
//         temp("i,j,c,d") = C("b,j") * corr("i,b,c,d");
//         corr("i,j,k,d") = C("c,k") * temp("i,j,c,d");
//         REQUIRE(tensorwrapper::tensor::allclose(rv, corr));
//     }

//     SECTION("Transform modes 0, 1, and 3") {
//         ao_map ao_spaces;
//         mo_map mo_spaces;
//         mo_spaces.emplace(0, mos);
//         mo_spaces.emplace(1, mos);
//         mo_spaces.emplace(3, mos);
//         ao_spaces.emplace(2, aos);
//         auto [rv] = mod.run_as<pt>(ao_spaces, mo_spaces, r12);
//         simde::type::tensor temp, corr;
//         corr("i,b,c,d") = C("a,i") * G("a,b,c,d");
//         temp("i,j,c,d") = C("b,j") * corr("i,b,c,d");
//         corr("i,j,c,k") = C("d,k") * temp("i,j,c,d");
//         REQUIRE(tensorwrapper::tensor::allclose(rv, corr));
//     }

//     SECTION("Transform modes 0, 2, and 3") {
//         ao_map ao_spaces;
//         mo_map mo_spaces;
//         mo_spaces.emplace(0, mos);
//         mo_spaces.emplace(2, mos);
//         mo_spaces.emplace(3, mos);
//         ao_spaces.emplace(1, aos);
//         auto [rv] = mod.run_as<pt>(ao_spaces, mo_spaces, r12);
//         simde::type::tensor temp, corr;
//         corr("i,b,c,d") = C("a,i") * G("a,b,c,d");
//         temp("i,b,j,d") = C("c,j") * corr("i,b,c,d");
//         corr("i,b,j,k") = C("d,k") * temp("i,b,j,d");
//         REQUIRE(tensorwrapper::tensor::allclose(rv, corr));
//     }

//     SECTION("Transform modes 1, 2, and 3") {
//         ao_map ao_spaces;
//         mo_map mo_spaces;
//         mo_spaces.emplace(1, mos);
//         mo_spaces.emplace(2, mos);
//         mo_spaces.emplace(3, mos);
//         ao_spaces.emplace(0, aos);
//         auto [rv] = mod.run_as<pt>(ao_spaces, mo_spaces, r12);
//         simde::type::tensor temp, corr;
//         corr("a,i,c,d") = C("b,i") * G("a,b,c,d");
//         temp("a,i,j,d") = C("c,j") * corr("a,i,c,d");
//         corr("a,i,j,k") = C("d,k") * temp("a,i,j,d");
//         REQUIRE(tensorwrapper::tensor::allclose(rv, corr));
//     }

//     SECTION("Transform all modes") {
//         ao_map ao_spaces;
//         mo_map mo_spaces;
//         mo_spaces.emplace(0, mos);
//         mo_spaces.emplace(1, mos);
//         mo_spaces.emplace(2, mos);
//         mo_spaces.emplace(3, mos);
//         auto [rv] = mod.run_as<pt>(ao_spaces, mo_spaces, r12);
//         simde::type::tensor corr, temp;
//         corr("i,b,c,d") = C("a,i") * G("a,b,c,d");
//         temp("i,j,c,d") = C("b,j") * corr("i,b,c,d");
//         corr("i,j,k,d") = C("c,k") * temp("i,j,c,d");
//         temp("i,j,k,l") = C("d,l") * corr("i,j,k,d");
//         REQUIRE(tensorwrapper::tensor::allclose(rv, temp));
//     }
// }

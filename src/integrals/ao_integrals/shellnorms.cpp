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

// #include "detail_/unpack_bases.hpp"
// #include "shellnorms.hpp"
// #include <simde/cauchy_schwarz_approximation.hpp>
// #include <simde/integral_factory.hpp>

// namespace integrals::ao_integrals {

// using factory_pt = simde::IntegralFactory;

// template<std::size_t NBodies, typename OperatorType>
// TEMPLATED_MODULE_CTOR(ShellNorms, NBodies, OperatorType) {
//     description("Calculates the Cauchy-Schwarz screening matrix for a pair of "
//                 "basis sets");

//     using my_pt = simde::ShellNorms<OperatorType>;
//     satisfies_property_type<my_pt>();

//     add_input<double>("Threshold")
//       .set_description("Convergence threshold of integrals")
//       .set_default(1.0E-16);

//     add_submodule<factory_pt>("AO Integral Factory")
//       .set_description("Used to generate the AO factory");
// }

// template<std::size_t NBodies, typename OperatorType>
// TEMPLATED_MODULE_RUN(ShellNorms, NBodies, OperatorType) {
//     using my_pt      = simde::ShellNorms<OperatorType>;
//     using elem_vec   = typename std::vector<double>;
//     using return_vec = typename std::vector<elem_vec>;

//     // Get inputs
//     auto bases     = detail_::unpack_bases<2>(inputs);
//     auto thresh    = inputs.at("Threshold").value<double>();
//     auto op_str    = OperatorType().as_string();
//     const auto& op = inputs.at(op_str).template value<const OperatorType&>();
//     const auto& fac_mod = submod.at("AO Integral Factory");

//     // Check if the basis sets are the same
//     bool same_bs = (bases[0] == bases[1]);

//     // Our return value
//     return_vec mat(bases[0].size(), elem_vec(bases[1].size(), 0.0));

//     // Lambda to fill in the values
//     std::function<void(int, int)> into_mat;
//     auto factory = fac_mod.run_as<factory_pt>(bases, op);
//     if constexpr(NBodies == 1) {
//         into_mat = [&, engine](int i, int j) mutable {
//             const auto& buf = factory.compute({i, j});
//             auto vals       = buf[0];

//             // Determine the number of compute values
//             std::size_t nvals = (bases[0][i].size() * bases[1][j].size());

//             // Find the norm and take the square root
//             double frobenius_norm_squared = 0.0;
//             if(vals != nullptr) {
//                 for(int a = 0; a < nvals; ++a) {
//                     frobenius_norm_squared += vals[a] * vals[a];
//                 }
//             }
//             mat[i][j] = std::sqrt(frobenius_norm_squared);
//             if(same_bs && (i != j)) {
//                 mat[j][i] = mat[i][j];
//             } // cut down on work
//         };
//     } else if constexpr(NBodies == 2) {
//         into_mat = [&, engine](int i, int j) mutable {
//             const auto& buf = = factory.compute({i, i, j, j});
//             auto vals         = buf[0];

//             // Determine the number of compute values
//             std::size_t nvals = (bases[0][i].size() * bases[0][i].size());
//             nvals *= (bases[1][j].size() * bases[1][j].size());

//             // Find the norm and take the square root
//             double inf_norm_squared = 0.0;
//             if(vals != nullptr) {
//                 for(int a = 0; a < nvals; ++a) {
//                     inf_norm_squared =
//                       std::max(inf_norm_squared, std::abs(vals[a]));
//                 }
//             }
//             mat[i][j] = std::sqrt(inf_norm_squared);
//             if(same_bs && (i != j)) {
//                 mat[j][i] = mat[i][j];
//             } // cut down on work
//         };
//     }

//     // Calculate values
//     auto& my_runtime = get_runtime();
//     auto& world      = my_runtime.madness_world();
//     for(int i = 0; i < bases[0].size(); ++i) {
//         // only do lower triangle if basis sets are the same
//         auto len = (same_bs) ? i : bases[1].size() - 1;
//         for(int j = 0; j <= len; ++j) { world.taskq.add(into_mat, i, j); }
//     }
//     world.gop.fence();

//     auto rv = results();
//     return my_pt::wrap_results(rv, mat);
// }

// template class ShellNorms<1, simde::type::el_identity>;
// template class ShellNorms<2, simde::type::el_el_coulomb>;
// template class ShellNorms<2, simde::type::el_el_stg>;
// template class ShellNorms<2, simde::type::el_el_yukawa>;

// } // namespace integrals::ao_integrals

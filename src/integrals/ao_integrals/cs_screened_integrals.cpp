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

// #include "cs_screened_integrals.hpp"
// #include "detail_/aos2shells.hpp"
// #include "detail_/get_coeff.hpp"
// #include "detail_/make_shape.hpp"
// #include "detail_/select_allocator.hpp"
// #include "detail_/shells2ord.hpp"
// #include "detail_/unpack_bases.hpp"
// #include <simde/cauchy_schwarz_approximation.hpp>
// #include <simde/integral_factory.hpp>
// #include <simde/tensor_representation/ao_tensor_representation.hpp>

// namespace integrals::ao_integrals {

// /// Grab the various detail_ functions
// using namespace detail_;

// using factory_pt = simde::IntegralFactory;

// /// For 1-body screening
// using identity_op_t = simde::type::el_identity;

// template<std::size_t N, typename OperatorType, bool direct>
// TEMPLATED_MODULE_CTOR(CSAOIntegral, N, OperatorType, direct) {
//     description("Computes an in-core integral with libint");
//     using my_pt = simde::AOTensorRepresentation<N, OperatorType>;

//     satisfies_property_type<my_pt>();

//     add_submodule<factory_pt>("AO Integral Factory")
//       .set_description("Used to generate the AO factory");

//     add_input<double>("Screening Threshold")
//       .set_default(0.0)
//       .set_description("Cauchy-Schwarz Screening Threshold");

//     if constexpr(N > 2) { /// Use operator based screening for 2-body integrals
//         add_submodule<simde::ShellNorms<OperatorType>>("Shell Norms")
//           .set_description(
//             "Computes the Cauchy-Schwarz Matrix for a pair of basis sets");
//     }

//     if constexpr(N == 2) { /// Use overlap based screening for 1-body integrals
//         add_submodule<simde::ShellNorms<identity_op_t>>("Shell Norms")
//           .set_description(
//             "Computes the Cauchy-Schwarz Matrix for a pair of basis sets");
//     }
// }

// template<std::size_t N, typename OperatorType, bool direct>
// TEMPLATED_MODULE_RUN(CSAOIntegral, N, OperatorType, direct) {
//     /// Typedefs
//     using my_pt         = simde::AOTensorRepresentation<N, OperatorType>;
//     using ao_space_t    = simde::type::ao_space;
//     using tensor_t      = simde::type::tensor;
//     using field_t       = typename tensor_t::field_type;
//     using size_vector_t = std::vector<std::size_t>;
//     using shell_norm_t  = std::vector<std::vector<double>>;

//     /// Grab input information
//     auto bases     = detail_::unpack_bases<N>(inputs);
//     auto op_str    = OperatorType().as_string();
//     auto cs_thresh = inputs.at("Screening Threshold").value<double>();
//     const auto& op = inputs.at(op_str).template value<const OperatorType&>();
//     const auto& fac_mod = submod.at("AO Integral Factory");

//     auto factory = fac_mod.run_as<factory_pt>(bases, op);
//     auto coeff   = detail_::get_coeff(op);

//     /// Calculate Shell Norms for screening
//     shell_norm_t mat1, mat2;
//     auto& cs_screen = submods.at("Shell Norms");

//     if constexpr(N == 2) {
//         auto bra = inputs.at("bra").template value<ao_space_t>();
//         auto ket = inputs.at("ket").template value<ao_space_t>();
//         identity_op_t I;
//         std::tie(mat1) =
//           cs_screen.run_as<simde::ShellNorms<identity_op_t>>(bra, I, ket);
//     }
//     if constexpr(N > 2) {
//         auto ket1 = inputs.at("ket 1").template value<ao_space_t>();
//         auto ket2 = inputs.at("ket 2").template value<ao_space_t>();
//         std::tie(mat1) =
//           cs_screen.run_as<simde::ShellNorms<OperatorType>>(ket1, op, ket2);
//     }
//     if constexpr(N == 4) {
//         auto bra1 = inputs.at("bra 1").template value<ao_space_t>();
//         auto bra2 = inputs.at("bra 2").template value<ao_space_t>();
//         std::tie(mat2) =
//           cs_screen.run_as<simde::ShellNorms<OperatorType>>(bra1, op, bra2);
//     }

//     /// Lambda to calculate values
//     auto l = [=](const auto& lo, const auto& up, auto* data) {
//         /// Convert index values from AOs to shells
//         size_vector_t lo_shells, up_shells;
//         for(auto i = 0; i < N; ++i) {
//             auto shells_in_tile = detail_::aos2shells(bases[i], lo[i], up[i]);
//             lo_shells.push_back(shells_in_tile.front());
//             up_shells.push_back(shells_in_tile.back());
//         }

//         /// Loop through shell combinations
//         size_vector_t curr_shells = lo_shells;
//         while(curr_shells[0] <= up_shells[0]) {
//             /// Determine which values will be computed this time
//             auto ord_pos =
//               detail_::shells2ord(bases, curr_shells, lo_shells, up_shells);

//             /// Determine the screening value
//             auto screen_value = mat1[curr_shells[N - 2]][curr_shells[N - 1]];
//             if constexpr(N == 4) {
//                 screen_value *= mat2[curr_shells[N - 4]][curr_shells[N - 3]];
//             }

//             /// Ensure that shells on the same center aren't screened out
//             /// in 1-body integrals
//             if constexpr(N == 2) {
//                 screen_value =
//                   (bases[0][curr_shells[0]].O == bases[1][curr_shells[1]].O) ?
//                     cs_thresh + 1 :
//                     screen_value;
//             }

//             /// If shell pair isn't screened out, calcuate the current values.
//             /// Otherwise, set the current values to 0.
//             if(screen_value > cs_thresh) {
//                 /// Compute values
//                 const auto& buf = factory.compute(curr_shells);
//                 auto vals       = buf[0];

//                 /// Copy libint values into tile data;
//                 for(auto i = 0; i < ord_pos.size(); ++i) {
//                     data[ord_pos[i]] = vals[i] * coeff;
//                 }
//             } else {
//                 for(auto i = 0; i < ord_pos.size(); ++i) {
//                     data[ord_pos[i]] = 0.0;
//                 }
//             }

//             /// Increment curr_shells
//             curr_shells[N - 1] += 1;
//             for(auto i = 1; i < N; ++i) {
//                 if(curr_shells[N - i] > up_shells[N - i]) {
//                     /// Reset this dimension and increment the next one
//                     /// curr_shells[0] accumulates until we reach the end
//                     curr_shells[N - i] = lo_shells[N - i];
//                     curr_shells[N - i - 1] += 1;
//                 }
//             }
//         }
//     };

//     tensor_t I(l, make_shape(bases),
//                select_allocator<direct, field_t>(bases, op, thresh, cs_thresh));

//     /// Finish
//     auto rv = results();
//     return my_pt::wrap_results(rv, I);
// }

// #define TEMPLATE_INT_AND_DIRECT(N, op)         \
//     template class CSAOintegral<N, op, false>; \
//     template class CSAOIntegral<N, op, true>

// TEMPLATE_INT_AND_DIRECT(2, simde::type::el_kinetic);
// TEMPLATE_INT_AND_DIRECT(2, simde::type::el_nuc_coulomb);
// TEMPLATE_INT_AND_DIRECT(2, simde::type::el_identity);
// TEMPLATE_INT_AND_DIRECT(3, simde::type::el_el_coulomb);
// TEMPLATE_INT_AND_DIRECT(4, simde::type::el_el_coulomb);
// TEMPLATE_INT_AND_DIRECT(3, simde::type::el_el_stg);
// TEMPLATE_INT_AND_DIRECT(4, simde::type::el_el_stg);
// TEMPLATE_INT_AND_DIRECT(3, simde::type::el_el_yukawa);
// TEMPLATE_INT_AND_DIRECT(4, simde::type::el_el_yukawa);

// #undef TEMPLATE_INT_AND_DIRECT

// } // namespace integrals::ao_integrals

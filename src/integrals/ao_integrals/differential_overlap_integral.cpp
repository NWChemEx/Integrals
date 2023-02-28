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

// #include "detail_/aos2shells.hpp"
// #include "detail_/get_coeff.hpp"
// #include "detail_/make_shape.hpp"
// #include "detail_/select_allocator.hpp"
// #include "detail_/shells2ord.hpp"
// #include "detail_/unpack_bases.hpp"
// #include <simde/tensor_representation/ao_tensor_representation.hpp>

// namespace integrals::ao_integrals {

// /// Grab the various detail_ functions
// using namespace detail_;

// using factory_pt = simde::IntegralFactory;

// template<bool direct>
// TEMPLATED_MODULE_CTOR(AOIntegralDOI, direct) {
//     description("Computes DOI integrals");
//     using my_pt = simde::AOTensorRepresentation<2, simde::type::el_el_delta>;

//     satisfies_property_type<my_pt>();

//     add_submodule<factory_pt>("AO Integral Factory")
//       .set_description("Used to generate the AO factory");
// }

// template<bool direct>
// TEMPLATED_MODULE_RUN(AOIntegralDOI, direct) {
//     using op_t          = simde::type::el_el_delta;
//     using my_pt         = simde::AOTensorRepresentation<2, op_t>;
//     using size_vector_t = std::vector<std::size_t>;
//     using tensor_t      = simde::type::tensor;
//     using field_t       = typename tensor_t::field_type;

//     auto init_bases     = detail_::unpack_bases<2>(inputs);
//     auto op_str         = op_t().as_string();
//     const auto& fac_mod = submod.at("AO Integral Factory");

//     auto factory = fac_mod.run_as<factory_pt>(bases, op);

//     /// Have to double up the basis sets
//     std::vector<simde::type::ao_basis_set> bases;
//     for(auto& set : init_bases) {
//         for(auto i = 0; i < 2; ++i) bases.push_back(set);
//     }

//     /// Lambda to calculate values
//     auto l = [=](const auto& lo, const auto& up, auto* data) {
//         /// Convert index values from AOs to shells
//         constexpr std::size_t N = 4;
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

//             /// Compute values
//             const auto& buf = factory.compute(curr_shells);
//             auto vals       = buf[0];

//             /// Copy libint values into tile data;
//             for(auto i = 0; i < ord_pos.size(); ++i) {
//                 data[ord_pos[i]] = vals[i];
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
//                select_allocator<direct, field_t>(bases, op, thresh));

//     /// Finish
//     auto rv = results();
//     return my_pt::wrap_results(rv, I);
// }

// template class AOIntegralDOI<false>;
// template class AOIntegralDOI<true>;

// } // namespace integrals::ao_integrals

// /*
//  * Copyright 2025 NWChemEx-Project
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

// #include "../utils/uncertainty_reductions.hpp"
// #include "ao_integrals.hpp"
// #include <integrals/integrals.hpp>
// #ifdef ENABLE_SIGMA
// #include <sigma/sigma.hpp>
// #endif

// using namespace tensorwrapper;

// namespace integrals::ao_integrals {
// namespace {

// template<typename FloatType, typename T, typename Tensor>
// auto average_error(T&& strides, T&& nbf, T&& ao_i, Tensor&& error,
//                    utils::mean_type mean) {
//     std::vector<FloatType> buffer;

//     for(std::size_t i = 0; i < nbf[0]; ++i) {
//         auto ioffset = (ao_i[0] + i) * strides[0];

//         for(std::size_t j = 0; j < nbf[1]; ++j) {
//             auto joffset = ioffset + (ao_i[1] + j) * strides[1];

//             for(std::size_t k = 0; k < nbf[2]; ++k) {
//                 auto koffset = joffset + (ao_i[2] + k) * strides[2];

//                 for(std::size_t l = 0; l < nbf[3]; ++l) {
//                     auto loffset = koffset + (ao_i[3] + l) * strides[3];
//                     buffer.push_back(error[loffset]);
//                 }
//             }
//         }
//     }

//     return utils::compute_mean(mean, buffer);
// }

// template<typename FloatType, typename T, typename Tensor>
// void update_block(T&& strides, T&& nbf, T&& ao_i, std::vector<FloatType>&
// out,
//                   Tensor&& value, FloatType error) {
//     for(std::size_t i = 0; i < nbf[0]; ++i) {
//         auto ioffset = (ao_i[0] + i) * strides[0];

//         for(std::size_t j = 0; j < nbf[1]; ++j) {
//             auto joffset = ioffset + (ao_i[1] + j) * strides[1];

//             for(std::size_t k = 0; k < nbf[2]; ++k) {
//                 auto koffset = joffset + (ao_i[2] + k) * strides[2];

//                 for(std::size_t l = 0; l < nbf[3]; ++l) {
//                     auto loffset = koffset + (ao_i[3] + l) * strides[3];
//                     out[loffset] = value[loffset] + error;
//                 }
//             }
//         }
//     }
// }

// template<typename FloatType, typename T, typename Tensor>
// void copy_block(T&& strides, T&& nbf, T&& ao_i, T&& new_nbf, T&& new_ao_i,
//                 std::vector<FloatType>& out) {
//     for(std::size_t i = 0; i < nbf[0]; ++i) {
//         auto ioffset = (ao_i[0] + i) * strides[0];

//         for(std::size_t j = 0; j < nbf[1]; ++j) {
//             auto joffset = ioffset + (ao_i[1] + j) * strides[1];

//             for(std::size_t k = 0; k < nbf[2]; ++k) {
//                 auto koffset = joffset + (ao_i[2] + k) * strides[2];

//                 for(std::size_t l = 0; l < nbf[3]; ++l) {
//                     auto loffset = koffset + (ao_i[3] + l) * strides[3];
//                     out[loffset] = value[loffset] + error;
//                 }
//             }
//         }
//     }
// }

// struct Kernel {
//     using shape_type = buffer::Contiguous::shape_type;
//     Kernel(shape_type shape, std::array<simde::type::ao_basis_set, 4> aos,
//            utils::mean_type mean) :
//       m_shape(std::move(shape)), m_aos(aos), m_mean(mean) {}

//     template<typename FloatType0, typename FloatType1>
//     Tensor operator()(const std::span<FloatType0> t,
//                       const std::span<FloatType1> error) {
//         throw std::runtime_error(
//           "UQ Integrals Driver kernel only supports same float types");
//     }

//     template<typename FloatType>
//     auto operator()(const std::span<FloatType> t,
//                     const std::span<FloatType> error) {
//         Tensor rv;

//         using float_type = std::decay_t<FloatType>;
//         if constexpr(types::is_uncertain_v<float_type>) {
//             throw std::runtime_error("Did not expect an uncertain type");
//         } else {
// #ifdef ENABLE_SIGMA
//             using tensorwrapper::buffer::make_contiguous;

//             std::array n_centers{m_aos[0].size(), m_aos[1].size(),
//                                  m_aos[2].size(), m_aos[3].size()};

//             std::array<std::size_t, 4> centers{0, 0, 0, 0};
//             std::array<std::size_t, 4> ao_i{0, 0, 0, 0};
//             std::array<std::size_t, 4> nbf{0, 0, 0, 0};

//             using uq_type = sigma::Uncertain<float_type>;
//             std::vector<uq_type> rv_data(m_shape.size());
//             std::array<std::size_t, 4> strides{0, 0, 0, 1};
//             strides[2] = strides[3] * m_aos[3].n_aos();
//             strides[1] = strides[2] * m_aos[2].n_aos();
//             strides[0] = strides[1] * m_aos[1].n_aos();

//             // If true, we can skip centers[1] > centers[0]
//             bool mu_is_nu = (m_aos[0] == m_aos[1]);
//             // If true, we can skip centers[3] > centers[2]
//             bool lam_is_sig = (m_aos[2] == m_aos[3]);
//             // If true, can skip centers[2] > centers[0] and centers[3] >
//             // centers[1]
//             bool all_same = (m_aos[0] == m_aos[2]) && mu_is_nu && lam_is_sig;

//             for(centers[0] = 0; centers[0] < n_centers[0]; ++centers[0]) {
//                 nbf[0] = m_aos[0][centers[0]].n_aos();

//                 ao_i[1] = 0;
//                 for(centers[1] = 0; centers[1] < n_centers[1]; ++centers[1])
//                 {
//                     if(centers[1] > centers[0] && mu_is_nu) break;
//                     nbf[1] = m_aos[1][centers[1]].n_aos();

//                     ao_i[2] = 0;
//                     for(centers[2] = 0; centers[2] < n_centers[2];
//                         ++centers[2]) {
//                         nbf[2] = m_aos[2][centers[2]].n_aos();

//                         ao_i[3] = 0;
//                         for(centers[3] = 0; centers[3] < n_centers[3];
//                             ++centers[3]) {
//                             if(centers[3] > centers[2] && lam_is_sig) break;

//                             nbf[3] = m_aos[3][centers[3]].n_aos();

//                             auto block_error = average_error<float_type>(
//                               strides, nbf, ao_i, error, m_mean);
//                             uq_type max_uq{0.0, block_error};

//                             update_block(strides, nbf, ao_i, rv_data, t,
//                                          max_uq);

//                             if(mu_is_nu) {
//                                 std::array<std::size_t, 4> new_nbf{
//                                   nbf[1], nbf[0], nbf[2], nbf[3]};
//                                 std::array<std::size_t, 4> new_ao_i{
//                                   ao_i[1], ao_i[0], ao_i[2], ao_i[3]};
//                                 copy_block(strides, nbf, ao_i, new_nbf,
//                                            new_ao_i, rv_data);
//                             }

//                             ao_i[3] += nbf[3];
//                         }
//                         ao_i[2] += nbf[2];
//                     }
//                     ao_i[1] += nbf[1];
//                 }
//                 ao_i[0] += nbf[0];
//             }
//             tensorwrapper::buffer::Contiguous t_w_contig(std::move(rv_data),
//                                                          m_shape);
//             rv = tensorwrapper::Tensor(m_shape, std::move(t_w_contig));
// #else
//             throw std::runtime_error("Sigma support not enabled!");
// #endif
//         }

//         return rv;
//     }
//     shape_type m_shape;
//     std::array<simde::type::ao_basis_set, 4> m_aos;
//     utils::mean_type m_mean;
// };

// const auto desc = R"(
// UQ Integrals Driver
// -------------------

// )";

// } // namespace

// using eri_pt   = simde::ERI4;
// using error_pt = integrals::property_types::Uncertainty<eri_pt>;

// MODULE_CTOR(UQAtomBlockedDriver) {
//     satisfies_property_type<eri_pt>();
//     description(desc);
//     add_submodule<eri_pt>("ERIs");
//     add_submodule<error_pt>("ERI Error");
//     add_input<std::string>("Mean Type").set_default("none");
// }

// MODULE_RUN(UQAtomBlockedDriver) {
//     const auto& [braket] = eri_pt::unwrap_inputs(inputs);
//     auto mean_str        = inputs.at("Mean Type").value<std::string>();
//     auto mean            = utils::mean_from_string(mean_str);

//     auto& eri_mod = submods.at("ERIs").value();
//     auto tol      = eri_mod.inputs().at("Threshold").value<double>();

//     const auto& t     = eri_mod.run_as<eri_pt>(braket);
//     const auto& error = submods.at("ERI Error").run_as<error_pt>(braket,
//     tol);

//     using tensorwrapper::buffer::make_contiguous;
//     const auto& t_buffer = make_contiguous(t.buffer());
//     const auto& e_buffer = make_contiguous(error.buffer());

//     const auto& bra = braket.bra();
//     const auto& ket = braket.ket();
//     const auto& mu  = bra.first.ao_basis_set();
//     const auto& nu  = bra.second.ao_basis_set();
//     const auto& lam = ket.first.ao_basis_set();
//     const auto& sig = ket.second.ao_basis_set();

//     std::array aos{mu, nu, lam, sig};

//     using buffer::visit_contiguous_buffer;
//     shape::Smooth shape =
//     t.buffer().layout().shape().as_smooth().make_smooth();

//     Kernel k(shape, aos, mean);
//     auto t_w_error = visit_contiguous_buffer(k, t_buffer, e_buffer);

//     auto rv = results();
//     return eri_pt::wrap_results(rv, t_w_error);
// }
// } // namespace integrals::ao_integrals

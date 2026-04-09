/*
 * Copyright 2025 NWChemEx-Project
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

#include "../utils/get_permutations.hpp"
#include "../utils/uncertainty_reductions.hpp"
#include "ao_integrals.hpp"
#include <cmath>
#include <integrals/integrals.hpp>
#ifdef ENABLE_SIGMA
#include <sigma/sigma.hpp>
#endif
using namespace tensorwrapper;

namespace integrals::ao_integrals {
namespace {

template<typename FloatType, typename T, typename Tensor>
auto average_error(T&& strides, T&& nbf, T&& ao_i, Tensor&& error,
                   utils::mean_type mean) {
#ifdef ENABLE_SIGMA
    using uq_type   = sigma::Interval<FloatType>;
    auto n_elements = nbf[0] * nbf[1] * nbf[2] * nbf[3];

    if(mean == utils::mean_type::none) {
        std::vector<uq_type> result;
        result.reserve(n_elements);
        for(std::size_t i = 0; i < nbf[0]; ++i) {
            auto ioffset = (ao_i[0] + i) * strides[0];
            for(std::size_t j = 0; j < nbf[1]; ++j) {
                auto joffset = ioffset + (ao_i[1] + j) * strides[1];
                for(std::size_t k = 0; k < nbf[2]; ++k) {
                    auto koffset = joffset + (ao_i[2] + k) * strides[2];
                    for(std::size_t l = 0; l < nbf[3]; ++l) {
                        auto loffset = koffset + (ao_i[3] + l) * strides[3];
                        FloatType w  = std::abs(error[loffset]);
                        result.push_back(uq_type{-w, w});
                    }
                }
            }
        }
        return result;
    }

    std::vector<FloatType> buffer;
    buffer.reserve(n_elements);
    for(std::size_t i = 0; i < nbf[0]; ++i) {
        auto ioffset = (ao_i[0] + i) * strides[0];
        for(std::size_t j = 0; j < nbf[1]; ++j) {
            auto joffset = ioffset + (ao_i[1] + j) * strides[1];
            for(std::size_t k = 0; k < nbf[2]; ++k) {
                auto koffset = joffset + (ao_i[2] + k) * strides[2];
                for(std::size_t l = 0; l < nbf[3]; ++l) {
                    auto loffset = koffset + (ao_i[3] + l) * strides[3];
                    buffer.push_back(error[loffset]);
                }
            }
        }
    }
    auto mean_value = utils::compute_mean(mean, buffer);
    FloatType w     = std::abs(mean_value);
    return std::vector<uq_type>(n_elements, uq_type{-w, w});
#else
    throw std::runtime_error("Sigma support not enabled!");
    return std::vector<int>{};
#endif
}

#ifdef ENABLE_SIGMA
template<typename FloatType, typename T, typename Tensor>
auto compute_block(T&& strides, T&& nbf, T&& ao_i, Tensor&& value,
                   const std::vector<sigma::Interval<FloatType>>& errors) {
    auto n_elements = nbf[0] * nbf[1] * nbf[2] * nbf[3];
    std::vector<sigma::Interval<FloatType>> buffer(n_elements);

    for(std::size_t i = 0; i < nbf[0]; ++i) {
        auto ilocal  = i * nbf[1] * nbf[2] * nbf[3];
        auto ioffset = (ao_i[0] + i) * strides[0];

        for(std::size_t j = 0; j < nbf[1]; ++j) {
            auto jlocal  = ilocal + j * nbf[2] * nbf[3];
            auto joffset = ioffset + (ao_i[1] + j) * strides[1];

            for(std::size_t k = 0; k < nbf[2]; ++k) {
                auto klocal  = jlocal + k * nbf[3];
                auto koffset = joffset + (ao_i[2] + k) * strides[2];

                for(std::size_t l = 0; l < nbf[3]; ++l) {
                    auto llocal    = klocal + l;
                    auto loffset   = koffset + (ao_i[3] + l) * strides[3];
                    buffer[llocal] = errors[llocal] + value[loffset];
                }
            }
        }
    }
    return buffer;
}
#endif

template<typename FloatType, typename T>
void set_block(T&& strides, T&& nbf,
               const std::array<std::size_t, 4>& permuted_ao_offsets,
               const std::array<std::size_t, 4>& sigma,
               const std::vector<FloatType>& block,
               std::vector<FloatType>& out) {
    // sigma[d] = the original mode that for what is now mode d. Therefore,
    // sigma[d] maps us back to the original mode, e.g., if the permutation
    // took 0, 1, 2, 3 to 3, 2, 1, 0 then sigma[0] = 3, sigma[1] = 2,
    // sigma[2] = 1, sigma[3] = 0.

    // If cidx is a 4-tuple of indices using the original modes, then for
    // output dimension d the new mode is cidx[sigma[d]].

    // Here we iterate in canonical (i,j,k,l) order — the same order the block
    // was filled — and then scatter to its permuted position in out.
    std::size_t block_idx = 0;
    for(std::size_t i = 0; i < nbf[0]; ++i) {
        for(std::size_t j = 0; j < nbf[1]; ++j) {
            for(std::size_t k = 0; k < nbf[2]; ++k) {
                for(std::size_t l = 0; l < nbf[3]; ++l) {
                    // This is the index using the original modes
                    std::array<std::size_t, 4> cidx{i, j, k, l};
                    const auto out_idx =
                      (permuted_ao_offsets[0] + cidx[sigma[0]]) * strides[0] +
                      (permuted_ao_offsets[1] + cidx[sigma[1]]) * strides[1] +
                      (permuted_ao_offsets[2] + cidx[sigma[2]]) * strides[2] +
                      (permuted_ao_offsets[3] + cidx[sigma[3]]) * strides[3];
                    out[out_idx] = block[block_idx++];
                }
            }
        }
    }
}

struct Kernel {
    using shape_type = buffer::Contiguous::shape_type;
    Kernel(shape_type shape, std::array<simde::type::ao_basis_set, 4> aos,
           utils::mean_type mean) :
      m_shape(std::move(shape)), m_aos(aos), m_mean(mean) {}

    template<typename FloatType0, typename FloatType1>
    Tensor operator()(const std::span<FloatType0> t,
                      const std::span<FloatType1> error) {
        throw std::runtime_error(
          "UQ Integrals Driver kernel only supports same float types");
    }

    template<typename FloatType>
    auto operator()(const std::span<FloatType> t,
                    const std::span<FloatType> error) {
        Tensor rv;

        using float_type = std::decay_t<FloatType>;
        if constexpr(types::is_uncertain_v<float_type> ||
                     types::is_interval_v<float_type>) {
            throw std::runtime_error("Did not expect an uncertain type");
        } else {
#ifdef ENABLE_SIGMA
            using tensorwrapper::buffer::make_contiguous;
            using utils::get_permutations_with_sigma;

            std::array n_centers{m_aos[0].size(), m_aos[1].size(),
                                 m_aos[2].size(), m_aos[3].size()};

            std::array<std::size_t, 4> centers{0, 0, 0, 0};
            std::array<std::size_t, 4> ao_offsets{0, 0, 0, 0};
            std::array<std::size_t, 4> nbf{0, 0, 0, 0};

            using uq_type = sigma::Interval<float_type>;
            std::vector<uq_type> rv_data(m_shape.size());
            std::array<std::size_t, 4> strides{0, 0, 0, 1};
            strides[2] = strides[3] * m_aos[3].n_aos();
            strides[1] = strides[2] * m_aos[2].n_aos();
            strides[0] = strides[1] * m_aos[1].n_aos();

            // If true, we can skip centers[1] > centers[0]
            bool mu_is_nu = (m_aos[0] == m_aos[1]);
            // If true, we can skip centers[3] > centers[2]
            bool lam_is_sig = (m_aos[2] == m_aos[3]);
            // If true, can skip centers[2] > centers[0] and centers[3] >
            // centers[1]
            bool all_same = (m_aos[0] == m_aos[2]) && mu_is_nu && lam_is_sig;

            for(centers[0] = 0; centers[0] < n_centers[0]; ++centers[0]) {
                nbf[0] = m_aos[0][centers[0]].n_aos();

                ao_offsets[1] = 0;
                for(centers[1] = 0; centers[1] < n_centers[1]; ++centers[1]) {
                    // We restrict our bra pairs to centers[0] <= centers[1]
                    if(centers[1] > centers[0] && mu_is_nu) break;
                    nbf[1] = m_aos[1][centers[1]].n_aos();

                    ao_offsets[2] = 0;
                    for(centers[2] = 0; centers[2] < n_centers[2];
                        ++centers[2]) {
                        // (c2, c3) <= (c0, c1) is impossible if c2 > c0
                        if(centers[2] > centers[0] && all_same) break;
                        bool c2eqc0 = centers[2] == centers[0];
                        nbf[2]      = m_aos[2][centers[2]].n_aos();

                        ao_offsets[3] = 0;
                        for(centers[3] = 0; centers[3] < n_centers[3];
                            ++centers[3]) {
                            // Restrict ket pairs to centers[2] <= centers[3]
                            if(centers[3] > centers[2] && lam_is_sig) break;

                            nbf[3] = m_aos[3][centers[3]].n_aos();
                            // Skip (c2,c3) > (c0,c1) lexicographically
                            bool pair_gt = (c2eqc0 && centers[3] > centers[1]);
                            if(pair_gt && all_same) break;

                            auto block_errors = average_error<float_type>(
                              strides, nbf, ao_offsets, error, m_mean);

                            // Compute (ab|cd)
                            auto block = compute_block(strides, nbf, ao_offsets,
                                                       t, block_errors);

                            // Set all symmetry equivalent blocks to `block`
                            auto perms = get_permutations_with_sigma(
                              ao_offsets, mu_is_nu, lam_is_sig, all_same);
                            for(auto& [perm, sigma] : perms) {
                                set_block(strides, nbf, perm, sigma, block,
                                          rv_data);
                            }

                            ao_offsets[3] += nbf[3];
                        }
                        ao_offsets[2] += nbf[2];
                    }
                    ao_offsets[1] += nbf[1];
                }
                ao_offsets[0] += nbf[0];
            }
            tensorwrapper::buffer::Contiguous t_w_contig(std::move(rv_data),
                                                         m_shape);
            rv = tensorwrapper::Tensor(m_shape, std::move(t_w_contig));
#else
            throw std::runtime_error("Sigma support not enabled!");
#endif
        }

        return rv;
    }
    shape_type m_shape;
    std::array<simde::type::ao_basis_set, 4> m_aos;
    utils::mean_type m_mean;
};

const auto desc = R"(
UQ Integrals Driver
-------------------

)";

} // namespace

using eri_pt   = simde::ERI4;
using error_pt = integrals::property_types::Uncertainty<eri_pt>;

MODULE_CTOR(UQAtomSymmBlockedDriver) {
    satisfies_property_type<eri_pt>();
    description(desc);
    add_submodule<eri_pt>("ERIs");
    add_submodule<error_pt>("ERI Error");
    add_input<std::string>("Mean Type").set_default("none");
}

MODULE_RUN(UQAtomSymmBlockedDriver) {
    const auto& [braket] = eri_pt::unwrap_inputs(inputs);
    auto mean_str        = inputs.at("Mean Type").value<std::string>();
    auto mean            = utils::mean_from_string(mean_str);

    auto& eri_mod = submods.at("ERIs").value();
    auto tol      = eri_mod.inputs().at("Threshold").value<double>();

    const auto& t     = eri_mod.run_as<eri_pt>(braket);
    const auto& error = submods.at("ERI Error").run_as<error_pt>(braket, tol);

    using tensorwrapper::buffer::make_contiguous;
    const auto& t_buffer = make_contiguous(t.buffer());
    const auto& e_buffer = make_contiguous(error.buffer());

    const auto& bra = braket.bra();
    const auto& ket = braket.ket();
    const auto& mu  = bra.first.ao_basis_set();
    const auto& nu  = bra.second.ao_basis_set();
    const auto& lam = ket.first.ao_basis_set();
    const auto& sig = ket.second.ao_basis_set();

    std::array aos{mu, nu, lam, sig};

    using buffer::visit_contiguous_buffer;
    shape::Smooth shape = t.buffer().layout().shape().as_smooth().make_smooth();

    Kernel k(shape, aos, mean);
    auto t_w_error = visit_contiguous_buffer(k, t_buffer, e_buffer);

    auto rv = results();
    return eri_pt::wrap_results(rv, t_w_error);
}
} // namespace integrals::ao_integrals

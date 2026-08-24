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
#include <integrals/integrals.hpp>
#ifdef ENABLE_SIGMA
#include <sigma/sigma.hpp>
#endif
using namespace tensorwrapper;

namespace integrals::ao_integrals {
namespace {

using wtf::buffer::FloatBuffer;

template<typename UQType, typename T, typename Tensor>
auto average_error(T&& strides, T&& nbf, T&& ao_i, Tensor&& error,
                   utils::mean_type mean,
                   const integrals::property_types::UQFactory& factory = {}) {
    std::string error_base =
      "integrals::ao_integrals::UQAtomSymmBlockedDriver: ";

#ifdef ENABLE_SIGMA
    using float_type = typename UQType::value_t;

    auto n_elements = nbf[0] * nbf[1] * nbf[2] * nbf[3];

    // The factory always returns a type-erased wtf::fp::Float (since it can
    // be configured, at runtime, to build any of the 5 UQ representations),
    // so every element is unwrapped to the concrete UQType requested by the
    // caller via float_cast.
    auto make_elem = [&](float_type ei) -> UQType {
        auto elem = factory(float_type(0.0), ei);
        return wtf::fp::float_cast<UQType>(elem);
    };

    if(mean == utils::mean_type::none) {
        FloatBuffer result;
        result.reserve<UQType>(n_elements);
        for(std::size_t i = 0; i < nbf[0]; ++i) {
            auto ioffset = (ao_i[0] + i) * strides[0];
            for(std::size_t j = 0; j < nbf[1]; ++j) {
                auto joffset = ioffset + (ao_i[1] + j) * strides[1];
                for(std::size_t k = 0; k < nbf[2]; ++k) {
                    auto koffset = joffset + (ao_i[2] + k) * strides[2];
                    for(std::size_t l = 0; l < nbf[3]; ++l) {
                        auto loffset = koffset + (ao_i[3] + l) * strides[3];
                        auto ei      = std::fabs(error[loffset]);
                        result.push_back<UQType>(make_elem(ei));
                    }
                }
            }
        }
        return result;
    }

    std::vector<float_type> buffer;
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
    auto value      = make_elem(mean_value);
    return FloatBuffer(std::vector<UQType>(n_elements, value));
#else
    throw std::runtime_error(error_base + "Sigma support not enabled!");
    return FloatBuffer{};
#endif
}

#ifdef ENABLE_SIGMA
template<typename UQType, typename T, typename Tensor>
auto compute_block(T&& strides, T&& nbf, T&& ao_i, Tensor&& value,
                   const FloatBuffer& errors) {
    auto n_elements  = nbf[0] * nbf[1] * nbf[2] * nbf[3];
    auto buffer      = FloatBuffer(std::vector<UQType>(n_elements));
    auto errors_span = errors.value<UQType>();
    auto buffer_span = buffer.value<UQType>();

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
                    auto llocal         = klocal + l;
                    auto loffset        = koffset + (ao_i[3] + l) * strides[3];
                    buffer_span[llocal] = errors_span[llocal] + value[loffset];
                }
            }
        }
    }
    return buffer;
}
#endif

template<typename UQType, typename T>
void set_block(T&& strides, T&& nbf,
               const std::array<std::size_t, 4>& permuted_ao_offsets,
               const std::array<std::size_t, 4>& sigma,
               const FloatBuffer& block, FloatBuffer& out) {
    // sigma[d] = the original mode that for what is now mode d. Therefore,
    // sigma[d] maps us back to the original mode, e.g., if the permutation
    // took 0, 1, 2, 3 to 3, 2, 1, 0 then sigma[0] = 3, sigma[1] = 2,
    // sigma[2] = 1, sigma[3] = 0.

    // If cidx is a 4-tuple of indices using the original modes, then for
    // output dimension d the new mode is cidx[sigma[d]].

    // Here we iterate in canonical (i,j,k,l) order — the same order the block
    // was filled — and then scatter to its permuted position in out.
    auto block_span       = block.value<UQType>();
    auto out_span         = out.value<UQType>();
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
                    out_span[out_idx] = block_span[block_idx++];
                }
            }
        }
    }
}

struct Kernel {
    using shape_type     = buffer::Contiguous::shape_type;
    using demangler_type = ::utilities::printing::Demangler;

    Kernel(shape_type shape, std::array<simde::type::ao_basis_set, 4> aos,
           utils::mean_type mean,
           integrals::property_types::UQFactory factory = {}) :
      m_shape(std::move(shape)),
      m_aos(aos),
      m_mean(mean),
      m_factory(std::move(factory)) {}

    template<typename FloatType0, typename FloatType1>
    Tensor operator()(const std::span<FloatType0> t,
                      const std::span<FloatType1> error) {
        auto type0 = demangler_type::demangle<FloatType0>();
        auto type1 = demangler_type::demangle<FloatType1>();
        throw std::runtime_error(
          m_error_base +
          "UQ Integrals Driver kernel only supports same "
          "float types. Got FloatType0 = " +
          type0 + " and FloatType1 = " + type1);
    }

    template<typename FloatType>
    auto operator()(const std::span<FloatType> t,
                    const std::span<FloatType> error) {
        Tensor rv;

        using float_type = std::decay_t<FloatType>;
        if constexpr(types::is_uq_type_v<float_type>) {
            auto type0 = demangler_type::demangle<FloatType>();
            throw std::runtime_error(
              m_error_base + "Expected non-UQ floating point type, but got " +
              type0);
        } else {
#ifdef ENABLE_SIGMA
            using utils::get_permutations_with_sigma;

            std::array n_centers{m_aos[0].size(), m_aos[1].size(),
                                 m_aos[2].size(), m_aos[3].size()};

            std::array<std::size_t, 4> centers{0, 0, 0, 0};
            std::array<std::size_t, 4> ao_offsets{0, 0, 0, 0};
            std::array<std::size_t, 4> nbf{0, 0, 0, 0};

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

            // The UQ kind is fixed by m_factory for the whole run (it never
            // varies block-to-block), so it's resolved to a concrete
            // UQType<float_type> exactly once here, via a single switch on
            // UQFactory::kind(). UQFactoryBase::operator() always returns a
            // wtf::fp::Float holding a *double*-valued UQ scalar (real ERI
            // tensors are always double in practice), so uq_type is built
            // from kind() + this Kernel's own float_type rather than by
            // inspecting a seed value, which would incorrectly decouple the
            // UQ scalar's width from float_type whenever this Kernel
            // template is instantiated for float.
            auto run_for_kind = [&]<template<typename> typename UQType> {
                using uq_type = UQType<float_type>;

                auto rv_data =
                  FloatBuffer(std::vector<uq_type>(m_shape.size()));

                for(centers[0] = 0; centers[0] < n_centers[0]; ++centers[0]) {
                    nbf[0] = m_aos[0][centers[0]].n_aos();

                    ao_offsets[1] = 0;
                    for(centers[1] = 0; centers[1] < n_centers[1];
                        ++centers[1]) {
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
                                // Restrict ket pairs to centers[2] <=
                                // centers[3]
                                if(centers[3] > centers[2] && lam_is_sig) break;

                                nbf[3] = m_aos[3][centers[3]].n_aos();
                                // Skip (c2,c3) > (c0,c1) lexicographically
                                bool pair_gt =
                                  (c2eqc0 && centers[3] > centers[1]);
                                if(pair_gt && all_same) break;

                                auto block_errors = average_error<uq_type>(
                                  strides, nbf, ao_offsets, error, m_mean,
                                  m_factory);

                                // Compute (ab|cd)
                                auto block = compute_block<uq_type>(
                                  strides, nbf, ao_offsets, t, block_errors);

                                // Set all symmetry equivalent blocks to
                                // `block`
                                auto perms = get_permutations_with_sigma(
                                  ao_offsets, mu_is_nu, lam_is_sig, all_same);
                                for(auto& [perm, sigma] : perms) {
                                    set_block<uq_type>(strides, nbf, perm,
                                                       sigma, block, rv_data);
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
            };

            using integrals::property_types::UQKind;
            switch(m_factory.kind()) {
                case UQKind::uncertain:
                    run_for_kind.template
                    operator()<tensorwrapper::types::uncertain_type>();
                    break;
                case UQKind::interval:
                    run_for_kind.template
                    operator()<tensorwrapper::types::interval_type>();
                    break;
                case UQKind::affine:
                    run_for_kind
                      .template operator()<tensorwrapper::types::affine_type>();
                    break;
                case UQKind::thresholded_affine:
                    run_for_kind.template
                    operator()<tensorwrapper::types::thresholded_affine_type>();
                    break;
                case UQKind::taylor_model:
                    run_for_kind.template
                    operator()<tensorwrapper::types::taylor_model_type>();
                    break;
            }
#else
            throw std::runtime_error(m_error_base +
                                     "Sigma support not enabled!");
#endif
        }

        return rv;
    }
    shape_type m_shape;
    std::array<simde::type::ao_basis_set, 4> m_aos;
    utils::mean_type m_mean;
    integrals::property_types::UQFactory m_factory;
    std::string m_error_base =
      "integrals::ao_integrals::UQAtomSymmBlockedDriver: ";
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
    add_submodule<integrals::property_types::UQInitializer>("UQ Initializer");
    add_input<std::string>("Mean Type").set_default("none");
}

MODULE_RUN(UQAtomSymmBlockedDriver) {
    const auto& [braket] = eri_pt::unwrap_inputs(inputs);
    auto mean_str        = inputs.at("Mean Type").value<std::string>();
    auto mean            = utils::mean_from_string(mean_str);

    auto factory = submods.at("UQ Initializer")
                     .run_as<integrals::property_types::UQInitializer>();

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

    Kernel k(shape, aos, mean, factory);
    simde::type::tensor t_w_error =
      visit_contiguous_buffer(k, t_buffer, e_buffer);

    auto rv = results();
    return eri_pt::wrap_results(rv, t_w_error);
}
} // namespace integrals::ao_integrals

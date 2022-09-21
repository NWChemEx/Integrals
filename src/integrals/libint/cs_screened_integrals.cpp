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

#include "cs_screened_integrals.hpp"
#include "detail_/aos2shells.hpp"
#include "detail_/bases_helper.hpp"
#include "detail_/make_engine.hpp"
#include "detail_/make_shape.hpp"
#include "detail_/select_allocator.hpp"
#include "detail_/shells2ord.hpp"
#include <simde/cauchy_schwarz_approximation.hpp>
#include <simde/tensor_representation/ao_tensor_representation.hpp>

namespace integrals {

/// Grab the various detail_ functions
using namespace detail_;

template<std::size_t N, typename OperatorType, bool direct>
TEMPLATED_MODULE_CTOR(CSLibint, N, OperatorType, direct) {
    description("Computes an in-core integral with libint");
    using my_pt = simde::AOTensorRepresentation<N, OperatorType>;

    satisfies_property_type<my_pt>();

    add_input<double>("Threshold")
      .set_default(1.0E-16)
      .set_description(
        "The target precision with which the integrals will be computed");

    add_input<double>("Screening Threshold")
      .set_default(0.0)
      .set_description("Cauchy-Schwarz Screening Threshold");

    using cs_approx_pt = simde::ShellNorms<OperatorType>;
    add_submodule<cs_approx_pt>("Shell Norms")
      .set_description(
        "Computes the Cauchy-Schwarz Matrix for a pair of basis sets");
}

template<std::size_t N, typename OperatorType, bool direct>
TEMPLATED_MODULE_RUN(CSLibint, N, OperatorType, direct) {
    /// Typedefs
    using my_pt         = simde::AOTensorRepresentation<N, OperatorType>;
    using cs_approx_pt  = simde::ShellNorms<OperatorType>;
    using size_vector_t = std::vector<std::size_t>;
    using tensor_t      = simde::type::tensor;
    using field_t       = typename tensor_t::field_type;
    using ao_space_t    = simde::type::ao_space;
    using shell_norm_t  = std::vector<std::vector<double>>;

    /// Grab input information
    auto bases     = unpack_bases<N>(inputs);
    auto op_str    = OperatorType().as_string();
    auto op        = inputs.at(op_str).template value<const OperatorType&>();
    auto thresh    = inputs.at("Threshold").value<double>();
    auto cs_thresh = inputs.at("Screening Threshold").value<double>();

    /// Calculate Shell Norms for screening
    shell_norm_t mat1, mat2;
    auto& cs_screen = submods.at("Shell Norms");
    if constexpr(N == 4) {
        auto bra1      = inputs.at("bra 1").template value<ao_space_t>();
        auto bra2      = inputs.at("bra 2").template value<ao_space_t>();
        std::tie(mat1) = cs_screen.run_as<cs_approx_pt>(bra1, op, bra2);
    }
    auto ket1      = inputs.at("ket 1").template value<ao_space_t>();
    auto ket2      = inputs.at("ket 2").template value<ao_space_t>();
    std::tie(mat2) = cs_screen.run_as<cs_approx_pt>(ket1, op, ket2);

    /// Geminal exponent handling
    constexpr auto is_stg =
      std::is_same_v<OperatorType, simde::type::el_el_stg>;
    constexpr auto is_yukawa =
      std::is_same_v<OperatorType, simde::type::el_el_yukawa>;

    double coeff = 1.0;
    if constexpr(is_stg || is_yukawa) {
        coeff = op.template at<0>().coefficient;
    }

    /// Lambda to calculate values
    auto l = [=](const auto& lo, const auto& up, auto* data) {
        /// Convert index values from AOs to shells
        size_vector_t lo_shells, up_shells;
        for(auto i = 0; i < N; ++i) {
            auto shells_in_tile = aos2shells(bases[i], lo[i], up[i]);
            lo_shells.push_back(shells_in_tile.front());
            up_shells.push_back(shells_in_tile.back());
        }

        /// Make the libint engine to calculate integrals
        auto engine     = make_engine(bases, op, thresh);
        const auto& buf = engine.results();

        /// Loop through shell combinations
        size_vector_t curr_shells = lo_shells;
        while(curr_shells[0] <= up_shells[0]) {
            /// Determine which values will be computed this time
            auto ord_pos = shells2ord(bases, curr_shells, lo_shells, up_shells);

            /// Check if current shells screen out
            auto screen_value = mat2[curr_shells[N - 2]][curr_shells[N - 1]];
            if constexpr(N == 4) {
                screen_value *= mat1[curr_shells[N - 4]][curr_shells[N - 3]];
            }

            if(screen_value > cs_thresh) {
                /// Compute values
                run_engine_(engine, bases, curr_shells,
                            std::make_index_sequence<N>());
                auto vals = buf[0];

                /// Copy libint values into tile data;
                for(auto i = 0; i < ord_pos.size(); ++i) {
                    data[ord_pos[i]] = vals[i] * coeff;
                }
            } else {
                for(auto i = 0; i < ord_pos.size(); ++i) {
                    data[ord_pos[i]] = 0.0;
                }
            }

            /// Increment curr_shells
            curr_shells[N - 1] += 1;
            for(auto i = 1; i < N; ++i) {
                if(curr_shells[N - i] > up_shells[N - i]) {
                    /// Reset this dimension and increment the next one
                    /// curr_shells[0] accumulates until we reach the end
                    curr_shells[N - i] = lo_shells[N - i];
                    curr_shells[N - i - 1] += 1;
                }
            }
        }
    };

    tensor_t I(l, make_shape(bases),
               select_allocator<direct, field_t>(bases, op, thresh, cs_thresh));

    /// Finish
    auto rv = results();
    return my_pt::wrap_results(rv, I);
}

#define TEMPLATE_INT_AND_DIRECT(N, op)     \
    template class CSLibint<N, op, false>; \
    template class CSLibint<N, op, true>

TEMPLATE_INT_AND_DIRECT(3, simde::type::el_el_coulomb);
TEMPLATE_INT_AND_DIRECT(4, simde::type::el_el_coulomb);
TEMPLATE_INT_AND_DIRECT(3, simde::type::el_el_stg);
TEMPLATE_INT_AND_DIRECT(4, simde::type::el_el_stg);
TEMPLATE_INT_AND_DIRECT(3, simde::type::el_el_yukawa);
TEMPLATE_INT_AND_DIRECT(4, simde::type::el_el_yukawa);

#undef TEMPLATE_INT_AND_DIRECT

} // namespace integrals

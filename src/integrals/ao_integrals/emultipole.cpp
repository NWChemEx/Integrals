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

#include "ao_integrals.hpp"
#include "detail_/aos2shells.hpp"
#include "detail_/bsets_shell_sizes.hpp"
#include "detail_/make_shape.hpp"
#include "detail_/select_allocator.hpp"
#include "detail_/shells2ord.hpp"
#include "detail_/unpack_bases.hpp"
#include <simde/integral_factory.hpp>
#include <simde/tensor_representation/ao_tensor_representation.hpp>

namespace integrals::ao_integrals {

using identity_op   = simde::type::el_identity;
using dipole_op     = simde::type::el_dipole;
using quadrupole_op = simde::type::el_quadrupole;
using octupole_op   = simde::type::el_octupole;

using overlap_pt    = simde::AOTensorRepresentation<2, identity_op>;
using dipole_pt     = simde::AOTensorRepresentation<2, dipole_op>;
using quadrupole_pt = simde::AOTensorRepresentation<2, quadrupole_op>;
using octupole_pt   = simde::AOTensorRepresentation<2, octupole_op>;

/// Type of a module that produces integral factories
template<typename OperatorType>
using factory_pt = simde::IntegralFactory<OperatorType>;
using factory_t  = simde::type::integral_factory;

/// Grab the various detail_ functions
using namespace detail_;

template<std::size_t L, typename OperatorType>
TEMPLATED_MODULE_CTOR(AOIntegralMultipole, L, OperatorType) {
    description("Computes an in-core multipole integral");

    /// This should satisfy overlap, but we can't reduce the dimensionality
    /// of the tensor at the moment.
    // satisfies_property_type<overlap_pt>();
    // identity_op I;
    // change_input(I.as_string()).change(std::move(I));

    satisfies_property_type<dipole_pt>();
    dipole_op r;
    change_input(r.as_string()).change(std::move(r));

    if constexpr(L > 0) {
        satisfies_property_type<quadrupole_pt>();
        quadrupole_op r2;
        change_input(r2.as_string()).change(std::move(r2));
    }

    if constexpr(L > 1) {
        satisfies_property_type<octupole_pt>();
        octupole_op r3;
        change_input(r3.as_string()).change(std::move(r3));
    }

    add_submodule<factory_pt<OperatorType>>("AO Integral Factory")
      .set_description("Used to generate the AO factory");
}

template<std::size_t L, typename OperatorType>
TEMPLATED_MODULE_RUN(AOIntegralMultipole, L, OperatorType) {
    // Typedefs
    using size_vector_t = std::vector<std::size_t>;
    using tensor_t      = simde::type::tensor;
    using field_t       = typename tensor_t::field_type;

    // Grab input information
    auto bases       = unpack_bases<2>(inputs);
    auto op_str      = OperatorType().as_string();
    auto op          = inputs.at(op_str).template value<const OperatorType&>();
    auto& fac_mod    = submods.at("AO Integral Factory");
    auto shell_sizes = bsets_shell_sizes(bases);

    // Cache result of factory module
    fac_mod.run_as<factory_pt<OperatorType>>(bases, op);

    // Lambda to calculate values
    auto l = [=](const auto& lo, const auto& up, auto* data) mutable {
        // Convert index values from AOs to shells
        // Leading index is for multipole components
        constexpr std::size_t N = 2;
        size_vector_t lo_shells, up_shells;
        for(auto i = 0; i < N; ++i) {
            auto shells_in_tile =
              aos2shells(shell_sizes[i], lo[i + 1], up[i + 1]);
            lo_shells.push_back(shells_in_tile.front());
            up_shells.push_back(shells_in_tile.back());
        }

        // Calculate the number of values per leading index
        auto leading_step = 0;
        for(auto i = lo_shells[0]; i <= up_shells[0]; ++i) {
            for(auto j = lo_shells[1]; j <= up_shells[1]; ++j) {
                leading_step += shell_sizes[0][i] * shell_sizes[1][j];
            }
        }

        // Get integral factory
        auto factory = fac_mod.run_as<factory_pt<OperatorType>>(bases, op);

        // Loop through shell combinations
        size_vector_t curr_shells = lo_shells;
        while(curr_shells[0] <= up_shells[0]) {
            // Determine which values will be computed this time
            auto ord_pos =
              shells2ord(shell_sizes, curr_shells, lo_shells, up_shells);

            // Compute values
            const auto& buf = factory.compute(curr_shells);

            // Copy libint values into tile data;
            for(auto i = lo[0]; i < up[0]; ++i) {
                auto depth = i * leading_step;
                for(auto j = 0; j < ord_pos.size(); ++j) {
                    data[ord_pos[j] + depth] = buf[i][j];
                }
            }

            // Increment curr_shells
            curr_shells[N - 1] += 1;
            for(auto i = 1; i < N; ++i) {
                if(curr_shells[N - i] > up_shells[N - i]) {
                    // Reset this dimension and increment the next one
                    // curr_shells[0] accumulates until we reach the end
                    curr_shells[N - i] = lo_shells[N - i];
                    curr_shells[N - i - 1] += 1;
                }
            }
        }
    };

    // Count up necessary components for multipole
    std::size_t leading_extent = 4;
    if constexpr(L == 1) { leading_extent = 10; }
    if constexpr(L == 2) { leading_extent = 20; }

    // Make complete tensor and slice out return values
    /// TODO: Switch out make_shape with an IntegralShape module
    tensor_t I(l, make_shape(bases, leading_extent),
               tensorwrapper::tensor::default_allocator<field_t>());

    auto rv = results();
    auto D  = I.slice({1, 0, 0}, {4, 7, 7});
    rv      = dipole_pt::wrap_results(rv, D);
    if constexpr(L > 0) {
        auto Q = I.slice({4, 0, 0}, {10, 7, 7});
        rv     = quadrupole_pt::wrap_results(rv, Q);
    }
    if constexpr(L > 1) {
        auto O = I.slice({10, 0, 0}, {20, 7, 7});
        rv     = octupole_pt::wrap_results(rv, O);
    }
    return rv;
}

// -----------------------------------------------------------------------------
// -- Template Declarations
// -----------------------------------------------------------------------------

template class AOIntegralMultipole<0, simde::type::el_dipole>;
template class AOIntegralMultipole<1, simde::type::el_quadrupole>;
template class AOIntegralMultipole<2, simde::type::el_octupole>;

} // namespace integrals::ao_integrals

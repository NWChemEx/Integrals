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

/// TODO: Optimize for IntegralFactory usage

#include "detail_/bsets_shell_sizes.hpp"
#include "detail_/unpack_bases.hpp"
#include "shellnorms.hpp"
#include <madness/world/MADworld.h>
#include <simde/cauchy_schwarz_approximation.hpp>
#include <simde/integral_factory.hpp>

namespace integrals::ao_integrals {

/// Type of a module that produces integral factories
template<typename OperatorType>
using factory_pt = simde::IntegralFactory<OperatorType>;
using factory_t  = simde::type::integral_factory;

/// Grab the various detail_ functions
using namespace detail_;

template<std::size_t NBodies, typename OperatorType>
TEMPLATED_MODULE_CTOR(ShellNorms, NBodies, OperatorType) {
    description("Calculates the Cauchy-Schwarz screening matrix for a pair of "
                "basis sets");

    using my_pt = simde::ShellNorms<OperatorType>;
    satisfies_property_type<my_pt>();

    add_submodule<factory_pt<OperatorType>>("AO Integral Factory")
      .set_description("Used to generate the AO factory");
}

template<std::size_t NBodies, typename OperatorType>
TEMPLATED_MODULE_RUN(ShellNorms, NBodies, OperatorType) {
    using my_pt      = simde::ShellNorms<OperatorType>;
    using elem_vec   = typename std::vector<double>;
    using return_vec = typename std::vector<elem_vec>;

    // Get inputs
    auto bases     = unpack_bases<2>(inputs);
    auto op_str    = OperatorType().as_string();
    auto& fac_mod  = submods.at("AO Integral Factory");
    const auto& op = inputs.at(op_str).template value<const OperatorType&>();

    // Double up on the basis sets if 2-body
    if constexpr(NBodies == 2) {
        bases.push_back(bases[0]);
        bases.push_back(bases[1]);
    }

    // Check if the basis sets are the same
    bool same_bs = (bases[0] == bases[1]);

    // Get shell sizes
    auto shell_sizes = bsets_shell_sizes(bases);

    // Our return value
    return_vec mat(bases[0].n_shells(), elem_vec(bases[1].n_shells(), 0.0));

    // Cache result of factory module
    fac_mod.run_as<factory_pt<OperatorType>>(bases, op);

    // Lambda to fill in the values
    std::function<void(std::size_t, std::size_t)> into_mat;
    if constexpr(NBodies == 1) {
        into_mat = [&](std::size_t i, std::size_t j) mutable {
            auto factory = fac_mod.run_as<factory_pt<OperatorType>>(bases, op);
            const auto& buf = factory.compute({i, j});
            auto vals       = buf[0];

            // Determine the number of compute values
            std::size_t nvals = (shell_sizes[0][i] * shell_sizes[1][j]);

            // Find the norm and take the square root
            double frobenius_norm_squared = 0.0;
            if(vals != nullptr) {
                for(int a = 0; a < nvals; ++a) {
                    frobenius_norm_squared += vals[a] * vals[a];
                }
            }
            mat[i][j] = std::sqrt(frobenius_norm_squared);
            if(same_bs && (i != j)) {
                mat[j][i] = mat[i][j];
            } // cut down on work
        };
    } else if constexpr(NBodies == 2) {
        into_mat = [&](std::size_t i, std::size_t j) mutable {
            auto factory = fac_mod.run_as<factory_pt<OperatorType>>(bases, op);
            const auto& buf = factory.compute({i, j, i, j});
            auto vals       = buf[0];

            // Determine the number of compute values
            std::size_t nvals = (shell_sizes[0][i] * shell_sizes[0][i]);
            nvals *= (shell_sizes[1][j] * shell_sizes[1][j]);

            // Find the norm and take the square root
            double inf_norm_squared = 0.0;
            if(vals != nullptr) {
                for(int a = 0; a < nvals; ++a) {
                    inf_norm_squared =
                      std::max(inf_norm_squared, std::abs(vals[a]));
                }
            }
            mat[i][j] = std::sqrt(inf_norm_squared);
            if(same_bs && (i != j)) {
                mat[j][i] = mat[i][j];
            } // cut down on work
        };
    }

    // Calculate values
    auto& my_runtime = get_runtime();
    auto comm        = get_runtime().mpi_comm();
    auto& world      = *madness::World::find_instance(SafeMPI::Intracomm(comm));
    for(std::size_t i = 0; i < bases[0].n_shells(); ++i) {
        // only do lower triangle if basis sets are the same
        auto len = (same_bs) ? i : bases[1].n_shells() - 1;
        for(std::size_t j = 0; j <= len; ++j) {
            world.taskq.add(into_mat, i, j);
        }
    }
    world.gop.fence();

    auto rv = results();
    return my_pt::wrap_results(rv, mat);
}

// -----------------------------------------------------------------------------
// -- Template Declarations
// -----------------------------------------------------------------------------

template class ShellNorms<1, simde::type::el_identity>;
template class ShellNorms<2, simde::type::el_el_coulomb>;
template class ShellNorms<2, simde::type::el_el_stg>;
template class ShellNorms<2, simde::type::el_el_yukawa>;

} // namespace integrals::ao_integrals

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
#include "cs_screened_integrals.hpp"
#include "detail_/aos2shells.hpp"
#include "detail_/get_coeff.hpp"
#include "detail_/make_shape.hpp"
#include "detail_/select_allocator.hpp"
#include "detail_/shells2ord.hpp"
#include "detail_/unpack_bases.hpp"
#include "shellnorms.hpp"
#include <simde/integral_factory.hpp>
#include <simde/tensor_representation/ao_tensor_representation.hpp>

namespace integrals::ao_integrals {

/// Type of a module that produces integral factories
template<typename OperatorType>
using factory_pt = simde::IntegralFactory<OperatorType>;
using factory_t  = simde::type::integral_factory;

/// Grab the various detail_ functions
using namespace detail_;

template<std::size_t N, typename OperatorType, bool direct>
TEMPLATED_MODULE_CTOR(AOIntegral, N, OperatorType, direct) {
    description("Computes integrals with Libint");
    using my_pt = simde::AOTensorRepresentation<N, OperatorType>;

    satisfies_property_type<my_pt>();

    add_submodule<factory_pt<OperatorType>>("AO Integral Factory")
      .set_description("Used to generate the AO factory");
}

template<std::size_t N, typename OperatorType, bool direct>
TEMPLATED_MODULE_RUN(AOIntegral, N, OperatorType, direct) {
    /// Typedefs
    using my_pt         = simde::AOTensorRepresentation<N, OperatorType>;
    using size_vector_t = std::vector<std::size_t>;
    using tensor_t      = simde::type::tensor;
    using field_t       = typename tensor_t::field_type;

    auto bases     = detail_::unpack_bases<N>(inputs);
    auto op_str    = OperatorType().as_string();
    auto& fac_mod  = submods.at("AO Integral Factory");
    const auto& op = inputs.at(op_str).template value<const OperatorType&>();

    auto [factory] = fac_mod.run_as<factory_pt<OperatorType>>(bases, op);
    auto coeff     = detail_::get_coefficient(op);

    /// Lambda to calculate values
    auto l = [=](const auto& lo, const auto& up, auto* data) mutable {
        /// Convert index values from AOs to shells
        size_vector_t lo_shells, up_shells;
        for(auto i = 0; i < N; ++i) {
            auto shells_in_tile = detail_::aos2shells(bases[i], lo[i], up[i]);
            lo_shells.push_back(shells_in_tile.front());
            up_shells.push_back(shells_in_tile.back());
        }

        /// Loop through shell combinations
        size_vector_t curr_shells = lo_shells;
        while(curr_shells[0] <= up_shells[0]) {
            /// Determine which values will be computed this time
            auto ord_pos =
              detail_::shells2ord(bases, curr_shells, lo_shells, up_shells);

            const auto& buf = factory.compute(curr_shells);
            auto vals       = buf[0];

            if(vals) {
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
               select_allocator<direct, field_t>(bases, op));

    /// Finish
    auto rv = results();
    return my_pt::wrap_results(rv, I);
}

// -----------------------------------------------------------------------------
// -- Template Declarations
// -----------------------------------------------------------------------------

#define TEMPLATE_INT_AND_DIRECT(N, op)       \
    template class AOIntegral<N, op, false>; \
    template class AOIntegral<N, op, true>

TEMPLATE_INT_AND_DIRECT(2, simde::type::el_el_coulomb);
TEMPLATE_INT_AND_DIRECT(3, simde::type::el_el_coulomb);
TEMPLATE_INT_AND_DIRECT(4, simde::type::el_el_coulomb);
TEMPLATE_INT_AND_DIRECT(2, simde::type::el_kinetic);
TEMPLATE_INT_AND_DIRECT(2, simde::type::el_nuc_coulomb);
TEMPLATE_INT_AND_DIRECT(2, simde::type::el_identity);
TEMPLATE_INT_AND_DIRECT(2, simde::type::el_el_stg);
TEMPLATE_INT_AND_DIRECT(3, simde::type::el_el_stg);
TEMPLATE_INT_AND_DIRECT(4, simde::type::el_el_stg);
TEMPLATE_INT_AND_DIRECT(2, simde::type::el_el_yukawa);
TEMPLATE_INT_AND_DIRECT(3, simde::type::el_el_yukawa);
TEMPLATE_INT_AND_DIRECT(4, simde::type::el_el_yukawa);
TEMPLATE_INT_AND_DIRECT(2, simde::type::el_el_f12_commutator);
TEMPLATE_INT_AND_DIRECT(3, simde::type::el_el_f12_commutator);
TEMPLATE_INT_AND_DIRECT(4, simde::type::el_el_f12_commutator);

#undef TEMPLATE_INT_AND_DIRECT

// -----------------------------------------------------------------------------
// -- Define Module Load Functions
// -----------------------------------------------------------------------------

#define ADD_INT_WITH_DIRECT(N, op, key_base)           \
    mm.add_module<AOIntegral<N, op, false>>(key_base); \
    mm.add_module<AOIntegral<N, op, true>>("Direct " key_base)

#define ADD_CS_INT_WITH_DIRECT(N, op, key_base)          \
    mm.add_module<CSAOIntegral<N, op, false>>(key_base); \
    mm.add_module<CSAOIntegral<N, op, true>>("Direct " key_base)

void load_ao_integrals(pluginplay::ModuleManager& mm) {
    // mm.add_module<AOIntegralDOI<false>>("DOI");
    // mm.add_module<AOIntegralDOI<true>>("Direct DOI");
    mm.add_module<AOIntegralMultipole<0, simde::type::el_dipole>>("EDipole");
    mm.add_module<AOIntegralMultipole<1, simde::type::el_quadrupole>>(
      "EQuadrupole");
    mm.add_module<AOIntegralMultipole<2, simde::type::el_octupole>>(
      "EOctupole");
    ADD_INT_WITH_DIRECT(2, simde::type::el_el_coulomb, "ERI2");
    ADD_INT_WITH_DIRECT(3, simde::type::el_el_coulomb, "ERI3");
    ADD_INT_WITH_DIRECT(4, simde::type::el_el_coulomb, "ERI4");
    ADD_INT_WITH_DIRECT(2, simde::type::el_kinetic, "Kinetic");
    ADD_INT_WITH_DIRECT(2, simde::type::el_nuc_coulomb, "Nuclear");
    ADD_INT_WITH_DIRECT(2, simde::type::el_identity, "Overlap");
    ADD_INT_WITH_DIRECT(2, simde::type::el_el_stg, "STG2");
    ADD_INT_WITH_DIRECT(3, simde::type::el_el_stg, "STG3");
    ADD_INT_WITH_DIRECT(4, simde::type::el_el_stg, "STG4");
    ADD_INT_WITH_DIRECT(2, simde::type::el_el_yukawa, "Yukawa2");
    ADD_INT_WITH_DIRECT(3, simde::type::el_el_yukawa, "Yukawa3");
    ADD_INT_WITH_DIRECT(4, simde::type::el_el_yukawa, "Yukawa4");
    // ADD_CS_INT_WITH_DIRECT(2, el_kinetic, "Kinetic CS");
    // ADD_CS_INT_WITH_DIRECT(2, el_nuc_coulomb, "Nuclear CS");
    // ADD_CS_INT_WITH_DIRECT(2, el_identity, "Overlap CS");
    // ADD_CS_INT_WITH_DIRECT(3, el_el_coulomb, "ERI3 CS");
    // ADD_CS_INT_WITH_DIRECT(4, el_el_coulomb, "ERI4 CS");
    // ADD_CS_INT_WITH_DIRECT(3, el_el_stg, "STG3 CS");
    // ADD_CS_INT_WITH_DIRECT(4, el_el_stg, "STG4 CS");
    // ADD_CS_INT_WITH_DIRECT(3, el_el_yukawa, "Yukawa3 CS");
    // ADD_CS_INT_WITH_DIRECT(4, el_el_yukawa, "Yukawa4 CS");

    // mm.add_module<ShellNormOverlap>("Shell Norms Overlap");
    // mm.add_module<ShellNormCoulomb>("Shell Norms Coulomb");
    // mm.add_module<ShellNormSTG>("Shell Norms STG");
    // mm.add_module<ShellNormYukawa>("Shell Norms Yukawa");

    // mm.add_module<AOIntegral<4, el_el_f12_commutator, false>>(
    //   "STG 4 Center dfdr Squared");
}

void ao_integrals_set_defaults(pluginplay::ModuleManager& mm) {
    // mm.change_submod("Kinetic CS", "Shell Norms", "Shell Norms Overlap");
    // mm.change_submod("Nuclear CS", "Shell Norms", "Shell Norms Overlap");
    // mm.change_submod("Overlap CS", "Shell Norms", "Shell Norms Overlap");
    // mm.change_submod("ERI3 CS", "Shell Norms", "Shell Norms Coulomb");
    // mm.change_submod("ERI4 CS", "Shell Norms", "Shell Norms Coulomb");
    // mm.change_submod("STG3 CS", "Shell Norms", "Shell Norms STG");
    // mm.change_submod("STG4 CS", "Shell Norms", "Shell Norms STG");
    // mm.change_submod("Yukawa3 CS", "Shell Norms", "Shell Norms Yukawa");
    // mm.change_submod("Yukawa4 CS", "Shell Norms", "Shell Norms Yukawa");
    // mm.change_submod("Direct Kinetic CS", "Shell Norms", "Shell Norms
    // Overlap"); mm.change_submod("Direct Nuclear CS", "Shell Norms", "Shell
    // Norms Overlap"); mm.change_submod("Direct Overlap CS", "Shell Norms",
    // "Shell Norms Overlap"); mm.change_submod("Direct ERI3 CS", "Shell Norms",
    // "Shell Norms Coulomb"); mm.change_submod("Direct ERI4 CS", "Shell Norms",
    // "Shell Norms Coulomb"); mm.change_submod("Direct STG3 CS", "Shell Norms",
    // "Shell Norms STG"); mm.change_submod("Direct STG4 CS", "Shell Norms",
    // "Shell Norms STG"); mm.change_submod("Direct Yukawa3 CS", "Shell Norms",
    // "Shell Norms Yukawa"); mm.change_submod("Direct Yukawa4 CS", "Shell
    // Norms", "Shell Norms Yukawa");
}

#undef ADD_INT_WITH_DIRECT
#undef ADD_CS_INT_WITH_DIRECT

} // namespace integrals::ao_integrals

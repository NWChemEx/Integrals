/*
 * Copyright 2026 NWChemEx-Project
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

#include "../utils/rank2_shell_norm.hpp"
#include "detail_/make_libint_basis_set.hpp"
#include "libint.hpp"
#include <cmath>
#include <integrals/property_types.hpp>

namespace integrals::libint {
namespace {

const auto desc = "";

} // namespace

using decontract_pt = integrals::property_types::DecontractBasisSet;
using eri4_pt       = simde::ERI4;
using pt            = integrals::property_types::PrimitivePairEstimator;

MODULE_CTOR(CauchySchwarzPrimitiveEstimator) {
    satisfies_property_type<pt>();
    description(desc);
    // TODO: Add citation for Chemist paper
    add_submodule<decontract_pt>("Decontract Basis Set");
    add_submodule<eri4_pt>("ERI4");
}

MODULE_RUN(CauchySchwarzPrimitiveEstimator) {
    const auto&& [bra_basis, ket_basis] = pt::unwrap_inputs(inputs);

    auto& to_prims_mod     = submods.at("Decontract Basis Set");
    const auto& bra_prims  = to_prims_mod.run_as<decontract_pt>(bra_basis);
    const auto& ket_prims  = to_prims_mod.run_as<decontract_pt>(ket_basis);
    const auto n_bra_prims = bra_prims.n_primitives();
    const auto n_ket_prims = ket_prims.n_primitives();

    // Should always be true, but we check for sanity
    assert(n_bra_prims == bra_basis.n_primitives());
    assert(n_ket_prims == ket_basis.n_primitives());

    // TODO: We only need the hyper diagonal, so this is very wasteful
    simde::type::aos_squared bra(bra_prims, ket_prims);
    simde::type::v_ee_type v_ee{};
    chemist::braket::BraKet mnls(bra, v_ee, bra);
    const auto& prim4 = submods.at("ERI4").run_as<eri4_pt>(mnls);

    // TODO: Make our basis set normalize itself.
    auto bra_libint = detail_::make_libint_basis_set(bra_basis);
    auto ket_libint = detail_::make_libint_basis_set(ket_basis);

    using tensorwrapper::buffer::make_contiguous;
    const auto& eris = make_contiguous(prim4.buffer());

    // TODO: Use floating point type of the basis sets
    using float_type = double;
    std::vector<float_type> data(n_bra_prims * n_ket_prims, 0.0);
    tensorwrapper::shape::Smooth shape({n_bra_prims, n_ket_prims});
    tensorwrapper::buffer::Contiguous buffer(std::move(data), shape);

    using iter_type    = std::decay_t<decltype(n_bra_prims)>; // Type of indices
    using index_array  = std::array<iter_type, 4>; // Type of a set of 4 indices
    using index_vector = std::vector<iter_type>; // Type of a vector of indices

    index_array ao_offsets{0, 0, 0, 0};
    index_array naos{0, 0, 0, 0};
    index_vector shell{0, 0};
    index_vector prim{0, 0};
    index_vector prim_offsets{0, 0};
    index_vector abs_prim{0, 0};

    for(shell[0] = 0; shell[0] < bra_basis.n_shells(); ++shell[0]) {
        const auto& bra_shell = bra_libint.at(shell[0]);
        assert(bra_shell.contr.size() == 1); // No general contraction support
        const auto& bra_coeff        = bra_shell.contr[0].coeff;
        const auto n_prims_bra_shell = bra_coeff.size();

        ao_offsets[0] = 0;
        ao_offsets[2] = 0;
        for(prim[0] = 0; prim[0] < n_prims_bra_shell; ++prim[0]) {
            const auto c_i = std::fabs(bra_coeff[prim[0]]);
            abs_prim[0]    = prim_offsets[0] + prim[0];
            naos[0]        = bra_basis.shell(shell[0]).size();
            naos[2]        = naos[0];

            prim_offsets[1] = 0;
            ao_offsets[1]   = 0;
            ao_offsets[3]   = 0;

            for(shell[1] = 0; shell[1] < ket_basis.n_shells(); ++shell[1]) {
                const auto& ket_shell = ket_libint.at(shell[1]);
                assert(ket_shell.contr.size() == 1); // No general contractions
                const auto& ket_coeff        = ket_shell.contr[0].coeff;
                const auto n_prims_ket_shell = ket_coeff.size();

                for(prim[1] = 0; prim[1] < n_prims_ket_shell; ++prim[1]) {
                    const auto c_j = std::fabs(ket_coeff[prim[1]]);
                    abs_prim[1]    = prim_offsets[1] + prim[1];

                    naos[1] = ket_basis.shell(shell[1]).size();
                    naos[3] = naos[1];

                    auto C_ij = c_i * c_j;

                    // ao_offsets/Naos needs to respectively be the offset for
                    // the first "AO" and the number of "AOs" in the
                    // decontracted ijij shell quartet

                    auto shell_norm =
                      utils::rank2_shell_norm(eris, ao_offsets, naos);
                    buffer.set_elem(abs_prim, C_ij * shell_norm);

                    ao_offsets[1] += naos[1];
                    ao_offsets[3] += naos[3];

                } // loop over ket primitives

                prim_offsets[1] += n_prims_ket_shell;
            } // loop over ket shells

            ao_offsets[0] += naos[0];
            ao_offsets[2] += naos[2];

        } // loop over bra primitives

        prim_offsets[0] += n_prims_bra_shell;

    } // loop over bra shells

    simde::type::tensor rv(shape, std::move(buffer));

    auto result = results();
    return pt::wrap_results(result, rv);
}

} // namespace integrals::libint

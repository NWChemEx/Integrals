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

#include "detail_/fill_tensor.hpp"
#include "detail_/get_basis_sets.hpp"
#include "libint.hpp"
#include <integrals/property_types.hpp>
#include <libint2.hpp>
namespace integrals::libint {
namespace {

const auto desc = R"(
Raw Primitive ERI4
==================

This module computes four-center electron repulsion integrals (ERIs) in a
fully decontracted (primitive) basis set with no libint normalization applied
to the shells.

The input braket's bra and ket basis sets are each decontracted via the
"Decontract Basis Set" submodule, which replaces every contracted shell with
one shell per primitive (coefficient set to 1.0). The resulting decontracted
basis sets are then passed to the "ERI4" submodule with libint's automatic
shell-coefficient normalization disabled, so the raw primitive integrals are
returned without any rescaling.

N.B. libint's normalization flag is a global static, so it is restored to
true immediately after the ERI4 call to avoid affecting other modules.
)";

} // namespace

using my_pt         = simde::ERI4;
using decontract_pt = integrals::property_types::DecontractBasisSet;

MODULE_CTOR(RawPrimitiveERIs) {
    satisfies_property_type<my_pt>();
    description(desc);
    add_submodule<decontract_pt>("Decontract Basis Set");
    add_input<double>("Threshold")
      .set_default(1.0E-16)
      .set_description(
        "The target precision with which the integrals will be computed");
    add_input<bool>("With UQ?").set_default(false);
}

MODULE_RUN(RawPrimitiveERIs) {
    const auto& [braket] = my_pt::unwrap_inputs(inputs);

    auto thresh  = inputs.at("Threshold").value<double>();
    auto with_uq = inputs.at("With UQ?").value<bool>();

    auto bra = braket.bra();
    auto ket = braket.ket();
    auto& op = braket.op();

    auto& decontract_mod = submods.at("Decontract Basis Set");

    // Decontract each of the four basis sets
    const auto& dc_bra0_bs =
      decontract_mod.run_as<decontract_pt>(bra.first.ao_basis_set());
    const auto& dc_bra1_bs =
      decontract_mod.run_as<decontract_pt>(bra.second.ao_basis_set());
    const auto& dc_ket0_bs =
      decontract_mod.run_as<decontract_pt>(ket.first.ao_basis_set());
    const auto& dc_ket1_bs =
      decontract_mod.run_as<decontract_pt>(ket.second.ao_basis_set());

    // Wrap the decontracted ao_basis_set objects back into aos/aos_squared
    simde::type::aos dc_bra0(dc_bra0_bs);
    simde::type::aos dc_bra1(dc_bra1_bs);
    simde::type::aos dc_ket0(dc_ket0_bs);
    simde::type::aos dc_ket1(dc_ket1_bs);

    simde::type::aos_squared dc_bra(dc_bra0, dc_bra1);
    simde::type::aos_squared dc_ket(dc_ket0, dc_ket1);

    // Disable libint's automatic shell normalization, compute, then restore
    auto original_normalization =
      libint2::Shell::do_enforce_unit_normalization();
    libint2::Shell::do_enforce_unit_normalization(false);
    auto& rv = this->get_runtime();

    auto basis_sets = detail_::get_basis_sets(dc_bra, dc_ket, false);
    auto t = detail_::fill_tensor<4>(basis_sets, op, rv, thresh, with_uq);
    libint2::Shell::do_enforce_unit_normalization(original_normalization);

    auto result = results();
    return my_pt::wrap_results(result, t);
}

} // namespace integrals::libint

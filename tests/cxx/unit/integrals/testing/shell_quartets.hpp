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

#pragma once
#include "ao_bases.hpp"
#include <simde/simde.hpp>

namespace integrals::testing {

/* After playing around, I found that if you set up a H2 dimer so that the
   H2 molecules are 3.0 Angstroms apart, then the (03|12) integral has a
   relatively large error (about 5.5E-9 using a tolerance of 1.0E-10). This
   function prepares 4, 1-shell AO basis sets so that you can easily evaluate
   that element.
 */
inline auto get_h2_dimer_0312_bases() {
    simde::type::ao_basis_set bra0, bra1, ket0, ket1;
    auto h2_2_aos = h2_sto3g_basis_set();

    // Centers 0 and 1 stay put, center 2 is 0 translated, and 3 is 1 translated
    bra0.add_center(h2_2_aos[0]); // 0
    bra1.add_center(h2_2_aos[1]); // 3
    ket0.add_center(h2_2_aos[1]); // 1
    ket1.add_center(h2_2_aos[0]); // 2

    const double r = 3.0;
    auto r0        = h2_2_aos[0].center();
    auto r1        = h2_2_aos[1].center();

    bra1[0].center().x() = r1.x() + r;
    bra1[0].center().y() = r1.y() + r;
    bra1[0].center().z() = r1.z() + r;
    ket1[0].center().x() = r0.x() + r;
    ket1[0].center().y() = r0.y() + r;
    ket1[0].center().z() = r0.z() + r;
    return std::make_tuple(bra0, bra1, ket0, ket1);
}

/* After playing around I found that my code was getting hard zero for this
 *  quartet, but libint got a non-zero, but small value.
 */
inline auto get_h2o_0034_bases() {
    simde::type::ao_basis_set bra0, ket3, ket4;
    auto water_aos = water_sto3g_basis_set();

    using abs_type    = simde::type::ao_basis_set::value_type;
    using cg_type     = simde::type::contracted_gaussian;
    using float_type  = double;
    using vector_type = std::vector<float_type>;

    const auto shell0   = water_aos.shell(0);
    const auto center0  = shell0.center().as_point();
    const auto cg0_view = shell0.contracted_gaussian();
    vector_type cs0{cg0_view[0].coefficient(), cg0_view[1].coefficient(),
                    cg0_view[2].coefficient()};
    vector_type es0{cg0_view[0].exponent(), cg0_view[1].exponent(),
                    cg0_view[2].exponent()};
    cg_type cg0(cs0.begin(), cs0.end(), es0.begin(), es0.end(), center0);

    abs_type abs0(center0);
    abs0.add_shell(shell0.pure(), shell0.l(), cg0);
    bra0.add_center(abs0);

    const auto shell3   = water_aos.shell(3);
    const auto center3  = shell3.center().as_point();
    const auto cg3_view = shell3.contracted_gaussian();
    vector_type cs3{cg3_view[0].coefficient(), cg3_view[1].coefficient(),
                    cg3_view[2].coefficient()};
    vector_type es3{cg3_view[0].exponent(), cg3_view[1].exponent(),
                    cg3_view[2].exponent()};
    cg_type cg3(cs3.begin(), cs3.end(), es3.begin(), es3.end(), center3);

    abs_type abs3(center3);
    abs3.add_shell(shell3.pure(), shell3.l(), cg3);
    ket3.add_center(abs3);

    const auto shell4   = water_aos.shell(4);
    const auto center4  = shell4.center().as_point();
    const auto cg4_view = shell4.contracted_gaussian();
    vector_type cs4{cg4_view[0].coefficient(), cg4_view[1].coefficient(),
                    cg4_view[2].coefficient()};
    vector_type es4{cg4_view[0].exponent(), cg4_view[1].exponent(),
                    cg4_view[2].exponent()};
    cg_type cg4(cs4.begin(), cs4.end(), es4.begin(), es4.end(), center4);

    abs_type abs4(center4);
    abs4.add_shell(shell4.pure(), shell4.l(), cg4);
    ket4.add_center(abs4);
    return std::make_tuple(bra0, bra0, ket3, ket4);
}

} // namespace integrals::testing

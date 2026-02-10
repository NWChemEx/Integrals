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

} // namespace integrals::testing
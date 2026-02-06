#pragma once
#include "testing.hpp"

namespace integrals::libint::test {

/* After playing around, I found that if you set up a H2 dimer so that the
   H2 molecules are 3.0 Angstroms apart, then the (03|12) integral has a
   relatively large error (about 5.5E-9 using a tolerance of 1.0E-10). This
   function prepares 4, 1-shell AO basis sets so that you can easily evaluate
   that element.
 */
inline auto get_h2_dimer_0312_bases() {
    simde::type::ao_basis_set bra0, bra1, ket0, ket1;
    auto h2_2_aos = ::test::h2_sto3g_basis_set();

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

/* Given four 1 shell basis sets, this function will create the shell quartet
   in the AO basis set by contracting the primitive integrals with the
   primitive coefficients. This is a sanity check more than anything else.
 */
template<typename BasisType>
inline auto manual_contract_shell(const BasisType& bra0, const BasisType& bra1,
                                  const BasisType& ket0, const BasisType& ket1,
                                  pluginplay::ModuleManager& mm) {
    auto dec_mod    = mm.at("Decontract Basis Set");
    using dec_pt    = integrals::property_types::DecontractBasisSet;
    auto bra0_prims = dec_mod.run_as<dec_pt>(bra0);
    auto bra1_prims = dec_mod.run_as<dec_pt>(bra1);
    auto ket0_prims = dec_mod.run_as<dec_pt>(ket0);
    auto ket1_prims = dec_mod.run_as<dec_pt>(ket1);
    simde::type::v_ee_type v_ee{};
    simde::type::aos_squared bra_prims(bra0_prims, bra1_prims);
    simde::type::aos_squared ket_prims(ket0_prims, ket1_prims);
    chemist::braket::BraKet mnls(bra_prims, v_ee, ket_prims);
    auto ints = mm.at("ERI4").run_as<simde::ERI4>(mnls);

    const auto& buffer = tensorwrapper::buffer::make_contiguous(ints.buffer());

    // TODO: expand to higher angular momentum
    double value = 0.0;
    for(std::size_t i = 0; i < bra0_prims.n_primitives(); ++i) {
        auto ci = bra0.primitive(i).coefficient();
        for(std::size_t j = 0; j < bra1_prims.n_primitives(); ++j) {
            auto cj = bra1.primitive(j).coefficient();
            for(std::size_t k = 0; k < ket0_prims.n_primitives(); ++k) {
                auto ck = ket0.primitive(k).coefficient();
                for(std::size_t l = 0; l < ket1_prims.n_primitives(); ++l) {
                    auto cl     = ket1.primitive(l).coefficient();
                    auto erased = buffer.get_elem({i, j, k, l});
                    auto val    = wtf::fp::float_cast<double>(erased);
                    value += ci * cj * ck * cl * val;
                }
            }
        }
    }
    return value;
}

} // namespace integrals::libint::test
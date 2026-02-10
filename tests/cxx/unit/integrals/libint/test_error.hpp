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

template<typename T>
double black_box_pair_metric(const T& prim_i, const T& prim_j) {
    auto ci       = prim_i.coefficient();
    auto cj       = prim_j.coefficient();
    auto zeta_i   = prim_i.exponent();
    auto zeta_j   = prim_j.exponent();
    auto dx       = prim_i.center().x() - prim_j.center().x();
    auto dy       = prim_i.center().y() - prim_j.center().y();
    auto dz       = prim_i.center().z() - prim_j.center().z();
    auto dr       = std::sqrt(dx * dx + dy * dy + dz * dz);
    auto exponent = (-zeta_i * zeta_j / (zeta_i + zeta_j)) * dr * dr;
    return ci * cj * std::exp(exponent);
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
    // auto K_mod    = mm.at("Black Box Primitive Pair Estimator");
    // using prim_pt = integrals::property_types::PrimitivePairEstimator;
    // auto K_ab     = K_mod.run_as<prim_pt>(bra0, bra1);
    // auto K_cd     = K_mod.run_as<prim_pt>(ket0, ket1);

    const auto& buffer = tensorwrapper::buffer::make_contiguous(ints.buffer());
    // const auto& K01    =
    // tensorwrapper::buffer::make_contiguous(K_ab.buffer()); const auto& K23 =
    // tensorwrapper::buffer::make_contiguous(K_cd.buffer());

    // TODO: expand to higher angular momentum
    double value = 0.0;
    double tol   = 1e-10;
    double error = 0.0;
    using wtf::fp::float_cast;
    for(std::size_t i = 0; i < bra0_prims.n_primitives(); ++i) {
        const auto prim_i = bra0.primitive(i);
        auto ci           = prim_i.coefficient();

        for(std::size_t j = 0; j < bra1_prims.n_primitives(); ++j) {
            const auto prim_j = bra1.primitive(j);
            auto cj           = prim_j.coefficient();

            auto kij = black_box_pair_metric(prim_i, prim_j);

            bool is01good = kij > tol;

            for(std::size_t k = 0; k < ket0_prims.n_primitives(); ++k) {
                const auto prim_k = ket0.primitive(k);
                auto ck           = prim_k.coefficient();

                for(std::size_t l = 0; l < ket1_prims.n_primitives(); ++l) {
                    const auto prim_l = ket1.primitive(l);
                    auto cl           = prim_l.coefficient();

                    auto kkl       = black_box_pair_metric(prim_k, prim_l);
                    bool is23good  = kkl > tol;
                    bool both_good = (kij * kkl) > tol;

                    auto erased = buffer.get_elem({i, j, k, l});
                    auto val    = wtf::fp::float_cast<double>(erased);
                    double Iijk = ci * cj * ck * cl * val;
                    if(is01good && is23good && both_good) {
                        value += Iijk;
                    } else {
                        error += std::fabs(Iijk);
                    }
                }
            }
        }
    }
    std::cout << "Computed error: " << simde::type::tensor(error) << std::endl;
    return simde::type::tensor(value);
}

} // namespace integrals::libint::test
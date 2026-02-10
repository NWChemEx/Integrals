#include "utils.hpp"

namespace integrals::utils {
namespace {

const auto desc = R"(
Primitive Contractor
====================

Given four AO basis sets, :math:`\{A\}, \{B\}, \{C\}, \{D\}`, this module will
build the ERIs in those basis sets. An element of the ERI,
:math:`(\mu\nu|\lambda\sigma)`, is built such that :math:`\mu \in \{A\}`,
:math:`\nu \in \{B\}`, :math:`\lambda \in \{C\}`, and :math:`\sigma \in \{D\}`.
)";

}

using decontract_pt = integrals::property_types::DecontractBasisSet;
using eris4_pt      = simde::ERI4;

MODULE_CTOR(PrimitiveContractor) {
    satisfies_property_type<eris4_pt>();
    description(desc);
    add_submodule<decontract_pt>("Decontract Basis Set")
      .set_description("Used to decontract the bases");
    add_submodule<eris4_pt>("Primitive ERI4")
      .set_description("Used to build the ERIs in the primitive bases");
}

MODULE_RUN(PrimitiveContractor) {
    const auto& [braket] = eris4_pt::unwrap_inputs(inputs);
    const auto& bra0     = braket.bra().first.ao_basis_set();
    const auto& bra1     = braket.bra().second.ao_basis_set();
    auto v_ee            = braket.op();
    const auto& ket0     = braket.ket().first.ao_basis_set();
    const auto& ket1     = braket.ket().second.ao_basis_set();

    // Step 1: Decontract the basis sets
    auto& dec_mod   = submods.at("Decontract Basis Set");
    auto bra0_prims = dec_mod.run_as<decontract_pt>(bra0);
    auto bra1_prims = dec_mod.run_as<decontract_pt>(bra1);
    auto ket0_prims = dec_mod.run_as<decontract_pt>(ket0);
    auto ket1_prims = dec_mod.run_as<decontract_pt>(ket1);

    // Step 2: Build the ERIs in the primitive basis sets
    simde::type::aos_squared bra_prims(bra0_prims, bra1_prims);
    simde::type::aos_squared ket_prims(ket0_prims, ket1_prims);
    chemist::braket::BraKet mnls(bra_prims, v_ee, ket_prims);
    auto ints = mm.at("Primitive ERI4").run_as<eris4_pt>(mnls);

    // TODO: Goes away when we can slice the tensor
    const auto& buffer = tensorwrapper::buffer::make_contiguous(ints.buffer());

    // // TODO: expand to higher angular momentum
    // double value = 0.0;
    // double tol   = 1e-10;
    // double error = 0.0;
    // using wtf::fp::float_cast;
    // for(std::size_t i = 0; i < bra0_prims.n_primitives(); ++i) {
    //     const auto prim_i = bra0.primitive(i);
    //     auto ci           = prim_i.coefficient();

    //     for(std::size_t j = 0; j < bra1_prims.n_primitives(); ++j) {
    //         const auto prim_j = bra1.primitive(j);
    //         auto cj           = prim_j.coefficient();

    //         auto kij = black_box_pair_metric(prim_i, prim_j);

    //         bool is01good = kij > tol;

    //         for(std::size_t k = 0; k < ket0_prims.n_primitives(); ++k) {
    //             const auto prim_k = ket0.primitive(k);
    //             auto ck           = prim_k.coefficient();

    //             for(std::size_t l = 0; l < ket1_prims.n_primitives(); ++l) {
    //                 const auto prim_l = ket1.primitive(l);
    //                 auto cl           = prim_l.coefficient();

    //                 auto kkl       = black_box_pair_metric(prim_k, prim_l);
    //                 bool is23good  = kkl > tol;
    //                 bool both_good = (kij * kkl) > tol;

    //                 auto erased = buffer.get_elem({i, j, k, l});
    //                 auto val    = wtf::fp::float_cast<double>(erased);
    //                 double Iijk = ci * cj * ck * cl * val;
    //                 if(is01good && is23good && both_good) {
    //                     value += Iijk;
    //                 } else {
    //                     error += std::fabs(Iijk);
    //                 }
    //             }
    //         }
    //     }
    // }
    // std::cout << "Computed error: " << simde::type::tensor(error) <<
    // std::endl; return simde::type::tensor(value);
}
} // namespace integrals::utils
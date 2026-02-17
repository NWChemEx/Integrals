#include "shell_quartet_iterator.hpp"
#include "utils.hpp"
#include <integrals/integrals.hpp>

namespace integrals::utils {
namespace {

/// Uses the black-box screening method to compute the AO element
template<typename ShellQuartet, typename OffsetType, typename TensorType>
decltype(auto) bb_ao_element(ShellQuartet&& shells, OffsetType&& prim_offsets,
                             OffsetType&& rel_ao, TensorType&& buffer,
                             double tol) {
    using float_type = double; // TODO: Get this from the buffer

    std::array<std::size_t, 4> prim{0, 0, 0, 0};
    std::vector<std::size_t> prim_index{0, 0, 0, 0};

    auto n_ao0   = shells[0].size();
    auto n_prim0 = shells[0].n_primitives();

    auto n_ao1   = shells[1].size();
    auto n_prim1 = shells[1].n_primitives();

    auto n_ao2   = shells[2].size();
    auto n_prim2 = shells[2].n_primitives();

    auto n_ao3    = shells[3].size();
    auto n_prims3 = shells[3].n_primitives();

    float_type error       = 0.0;
    float_type value       = 0.0;
    float_type bb_estimate = 0.0;

    for(prim[0] = 0; prim[0] < n_prim0; ++prim[0]) {
        auto c0       = shells[0].primitive(prim[0]).coefficient();
        auto z0       = shells[0].primitive(prim[0]).exponent();
        auto r0       = shells[0].primitive(prim[0]).center();
        prim_index[0] = prim_offsets[0] + prim[0] * n_ao0 + rel_ao[0];

        for(prim[1] = 0; prim[1] < n_prim1; ++prim[1]) {
            auto c1       = shells[1].primitive(prim[1]).coefficient();
            auto z1       = shells[1].primitive(prim[1]).exponent();
            auto r1       = shells[1].primitive(prim[1]).center();
            prim_index[1] = prim_offsets[1] + prim[1] * n_ao1 + rel_ao[1];

            auto dr01         = (r0 - r1).magnitude();
            auto z01          = -(z0 * z1) / (z0 + z1);
            auto k01          = c0 * c1 * std::exp(z01 * dr01 * dr01);
            bool k01_is_small = std::abs(k01) < tol;

            for(prim[2] = 0; prim[2] < n_prim2; ++prim[2]) {
                auto c2       = shells[2].primitive(prim[2]).coefficient();
                auto z2       = shells[2].primitive(prim[2]).exponent();
                auto r2       = shells[2].primitive(prim[2]).center();
                prim_index[2] = prim_offsets[2] + prim[2] * n_ao2 + rel_ao[2];

                for(prim[3] = 0; prim[3] < n_prims3; ++prim[3]) {
                    auto c3 = shells[3].primitive(prim[3]).coefficient();
                    auto z3 = shells[3].primitive(prim[3]).exponent();
                    auto r3 = shells[3].primitive(prim[3]).center();
                    prim_index[3] =
                      prim_offsets[3] + prim[3] * n_ao3 + rel_ao[3];

                    auto dr23          = (r2 - r3).magnitude();
                    auto z23           = -(z2 * z3) / (z2 + z3);
                    auto k23           = c2 * c3 * std::exp(z23 * dr23 * dr23);
                    bool k23_is_small  = std::abs(k23) < tol;
                    bool prod_is_small = std::abs(k01 * k23) < tol;

                    auto erased_elem = buffer.get_elem(prim_index);
                    auto elem = wtf::fp::float_cast<float_type>(erased_elem);

                    if(k01_is_small || k23_is_small || prod_is_small) {
                        error += c0 * c1 * c2 * c3 * elem;
                        bb_estimate += std::abs(k01 * k23);
                    } else {
                        value += c0 * c1 * c2 * c3 * elem;
                    }
                }
            }
        }
    }
    return std::make_tuple(value, error, bb_estimate);
}

/// Computes all AO quartets in a shell quartet and puts them in `result`.
template<typename ShellQuartet, typename OffsetType, typename TensorType1,
         typename TensorType2>
void fill_ao_quartet(ShellQuartet&& shells, OffsetType&& ao_offsets,
                     OffsetType&& prim_offsets, TensorType1&& result,
                     TensorType2&& buffer, TensorType1&& error_tensor,
                     TensorType1&& bb_estimate_tensor, double tol) {
    auto n_ao0 = shells[0].size();
    auto n_ao1 = shells[1].size();
    auto n_ao2 = shells[2].size();
    auto n_ao3 = shells[3].size();

    std::array<std::size_t, 4> ao{0, 0, 0, 0};
    std::vector<std::size_t> ao_index{0, 0, 0, 0};

    for(ao[0] = 0; ao[0] < n_ao0; ++ao[0]) {
        ao_index[0] = ao_offsets[0] + ao[0];

        for(ao[1] = 0; ao[1] < n_ao1; ++ao[1]) {
            ao_index[1] = ao_offsets[1] + ao[1];

            for(ao[2] = 0; ao[2] < n_ao2; ++ao[2]) {
                ao_index[2] = ao_offsets[2] + ao[2];

                for(ao[3] = 0; ao[3] < n_ao3; ++ao[3]) {
                    ao_index[3] = ao_offsets[3] + ao[3];

                    auto x =
                      bb_ao_element(shells, prim_offsets, ao, buffer, tol);
                    result.set_elem(ao_index, std::get<0>(x));
                    error_tensor.set_elem(ao_index, std::get<1>(x));
                    bb_estimate_tensor.set_elem(ao_index, std::get<2>(x));
                }
            }
        }
    }
}

const auto desc = R"(
Primitive Contractor
====================

Given four AO basis sets, :math:`\{A\}, \{B\}, \{C\}, \{D\}`, this module will
build the ERIs in those basis sets. An element of the ERI,
:math:`(\mu\nu|\lambda\sigma)`, is built such that :math:`\mu \in \{A\}`,
:math:`\nu \in \{B\}`, :math:`\lambda \in \{C\}`, and :math:`\sigma \in \{D\}`.
)";
} // namespace

using decontract_pt = integrals::property_types::DecontractBasisSet;
using eris4_pt      = simde::ERI4;

MODULE_CTOR(PrimitiveContractor) {
    satisfies_property_type<eris4_pt>();
    description(desc);
    add_submodule<decontract_pt>("Decontracter")
      .set_description("Used to decontract the bases");
    add_submodule<eris4_pt>("Primitive ERI4")
      .set_description("Used to build the ERIs in the primitive bases");
    add_input<double>("threshold");
    add_result<simde::type::tensor>("Corr ERI4");
    add_result<simde::type::tensor>("Actual Error");
    add_result<simde::type::tensor>("Black-Box Error Estimate");
}

MODULE_RUN(PrimitiveContractor) {
    using float_type = double; // TODO: Get this from the buffer

    const auto& [braket] = eris4_pt::unwrap_inputs(inputs);
    const auto tol       = inputs.at("threshold").value<double>();
    const auto& bra0     = braket.bra().first.ao_basis_set();
    const auto& bra1     = braket.bra().second.ao_basis_set();
    auto v_ee            = braket.op();
    const auto& ket0     = braket.ket().first.ao_basis_set();
    const auto& ket1     = braket.ket().second.ao_basis_set();

    // Step 1: Decontract the basis sets
    auto& dec_mod   = submods.at("Decontracter");
    auto bra0_prims = dec_mod.run_as<decontract_pt>(bra0);
    auto bra1_prims = dec_mod.run_as<decontract_pt>(bra1);
    auto ket0_prims = dec_mod.run_as<decontract_pt>(ket0);
    auto ket1_prims = dec_mod.run_as<decontract_pt>(ket1);

    // Step 2: Build the ERIs in the primitive basis sets
    simde::type::aos_squared bra_prims(bra0_prims, bra1_prims);
    simde::type::aos_squared ket_prims(ket0_prims, ket1_prims);
    chemist::braket::BraKet mnls(bra_prims, v_ee, ket_prims);
    auto ints = submods.at("Primitive ERI4").run_as<eris4_pt>(mnls);

    // TODO: Goes away when we can slice the tensor
    const auto& buffer = tensorwrapper::buffer::make_contiguous(ints.buffer());

    // Step 3: Make the result tensor
    tensorwrapper::shape::Smooth shape{bra0.n_aos(), bra1.n_aos(), ket0.n_aos(),
                                       ket1.n_aos()};
    std::vector<float_type> data(shape.size());
    tensorwrapper::buffer::Contiguous result_buffer(std::move(data), shape);
    simde::type::tensor tw_result(shape, std::move(result_buffer));
    auto& result = tensorwrapper::buffer::make_contiguous(tw_result.buffer());

    simde::type::tensor tw_error(tw_result);
    auto& error_buffer =
      tensorwrapper::buffer::make_contiguous(tw_error.buffer());

    simde::type::tensor tw_bb_estimate(tw_result);
    auto& bb_estimate_buffer =
      tensorwrapper::buffer::make_contiguous(tw_bb_estimate.buffer());

    // Step 4: Quantities needed for the for loops.
    using offset_type = std::array<std::size_t, 4>;
    offset_type shell{0, 0, 0, 0};
    offset_type nshells{bra0.n_shells(), bra1.n_shells(), ket0.n_shells(),
                        ket1.n_shells()};

    offset_type ao_offsets{0, 0, 0, 0};
    offset_type prim_offsets{0, 0, 0, 0};

    using shell_view = std::decay_t<decltype(bra0.shell(0))>;
    std::array<shell_view, 4> shell_views;

    // Step 5: Loop over the shell quartets and fill in the AO quartets
    for(shell[0] = 0; shell[0] < nshells[0]; ++shell[0]) {
        shell_views[0] = bra0.shell(shell[0]);
        auto n_ao0     = shell_views[0].size();
        auto n_prim0   = shell_views[0].n_primitives();

        ao_offsets[1]   = 0;
        prim_offsets[1] = 0;
        for(shell[1] = 0; shell[1] < nshells[1]; ++shell[1]) {
            shell_views[1] = bra1.shell(shell[1]);
            auto n_ao1     = shell_views[1].size();
            auto n_prim1   = shell_views[1].n_primitives();

            ao_offsets[2]   = 0;
            prim_offsets[2] = 0;
            for(shell[2] = 0; shell[2] < nshells[2]; ++shell[2]) {
                shell_views[2] = ket0.shell(shell[2]);
                auto n_ao2     = shell_views[2].size();
                auto n_prim2   = shell_views[2].n_primitives();

                ao_offsets[3]   = 0;
                prim_offsets[3] = 0;
                for(shell[3] = 0; shell[3] < nshells[3]; ++shell[3]) {
                    shell_views[3] = ket1.shell(shell[3]);
                    auto n_ao3     = shell_views[3].size();
                    auto n_prims3  = shell_views[3].n_primitives();

                    fill_ao_quartet(shell_views, ao_offsets, prim_offsets,
                                    result, buffer, error_buffer,
                                    bb_estimate_buffer, tol);

                    ao_offsets[3] += n_ao3;
                    prim_offsets[3] += n_prims3 * n_ao3;
                }

                ao_offsets[2] += n_ao2;
                prim_offsets[2] += n_prim2 * n_ao2;
            }

            ao_offsets[1] += n_ao1;
            prim_offsets[1] += n_prim1 * n_ao1;
        }
        ao_offsets[0] += n_ao0;
        prim_offsets[0] += n_prim0 * n_ao0;
    }

    simde::type::tensor corr;
    corr("i,j,k,l") = tw_result("i,j,k,l") + tw_error("i,j,k,l");

    auto rv = results();
    rv["Corr ERI4"].change(corr);
    rv["Actual Error"].change(tw_error);
    rv["Black-Box Error Estimate"].change(tw_bb_estimate);
    return eris4_pt::wrap_results(rv, tw_result);
}
} // namespace integrals::utils
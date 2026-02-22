#include "libint.hpp"
#include <array>
#include <integrals/integrals.hpp>

namespace integrals::libint {
namespace {

const auto desc = "Compute uncertainty estimates for ERI4 integrals";

// Sums contributions Q_abQ_cd when Q_ab or Q_cd is below tol
template<typename ShellsType, typename OffsetTypes, typename TensorType>
auto shell_block_error(ShellsType&& shells, OffsetTypes&& prim_offsets,
                       TensorType&& Q_ab, TensorType&& Q_cd, TensorType&& K_ab,
                       TensorType&& K_cd, double tol) {
    // Get the number of primitives in each shell
    const auto n_prim0 = shells[0].n_primitives();
    const auto n_prim1 = shells[1].n_primitives();
    const auto n_prim2 = shells[2].n_primitives();
    const auto n_prim3 = shells[3].n_primitives();

    using float_type = double; // TODO: get from Q_ab or Q_cd
    using iter_type = std::decay_t<decltype(n_prim0)>; // Type of the loop index

    // Create an array to hold the primitive indices
    std::array<iter_type, 4> prim{0, 0, 0, 0};

    // This is the accumulated error for this shell block
    double error = 0.0;

    using wtf::fp::float_cast;
    for(prim[0] = 0; prim[0] < n_prim0; ++prim[0]) {
        const auto p0 = prim_offsets[0] + prim[0];

        for(prim[1] = 0; prim[1] < n_prim1; ++prim[1]) {
            const auto p1 = prim_offsets[1] + prim[1];

            // Was the Q_ab element neglected for primitive pair (p0, p1)?
            std::vector i01{p0, p1};

            auto K_01           = float_cast<float_type>(K_ab.get_elem(i01));
            auto K_01_mag       = std::fabs(K_01);
            bool K_01_neglected = K_01_mag < tol;
            auto Q_01           = float_cast<float_type>(Q_ab.get_elem(i01));
            auto Q_01_mag       = std::fabs(Q_01);

            for(prim[2] = 0; prim[2] < n_prim2; ++prim[2]) {
                const auto p2 = prim_offsets[2] + prim[2];

                for(prim[3] = 0; prim[3] < n_prim3; ++prim[3]) {
                    const auto p3 = prim_offsets[3] + prim[3];

                    // Was the Q_cd element neglected for pair (p2, p3)?
                    std::vector i23{p2, p3};
                    auto K_23     = float_cast<float_type>(K_cd.get_elem(i23));
                    auto K_23_mag = std::fabs(K_23);
                    bool K_23_neglected = K_23_mag < tol;
                    auto Q_23     = float_cast<float_type>(Q_cd.get_elem(i23));
                    auto Q_23_mag = std::fabs(Q_23);

                    auto prod_neglected = K_01_mag * K_23_mag < tol;

                    // If either was neglected we pick up an error Q_01 * Q_23
                    if(K_01_neglected || K_23_neglected || prod_neglected) {
                        error += Q_01_mag * Q_23_mag;
                    }
                } // End loop over prim[3]
            } // End loop over prim[2]
        } // End loop over prim[1]
    } // End loop over prim[0]

    return error;
}

// Sets all AO integrals in the shell block to the computed error
template<typename ShellsType, typename OffsetTypes, typename TensorType>
void fill_ao_block(ShellsType&& shells, OffsetTypes&& ao_offsets,
                   TensorType&& t, double error) {
    // Get the number of AOs in each shell
    const auto n_aos0 = shells[0].size();
    const auto n_aos1 = shells[1].size();
    const auto n_aos2 = shells[2].size();
    const auto n_aos3 = shells[3].size();

    // Make an array to hold the AO indices
    using iter_type = std::decay_t<decltype(n_aos0)>; // Type of the loop index
    std::array<iter_type, 4> ao{0, 0, 0, 0};

    for(ao[0] = 0; ao[0] < n_aos0; ++ao[0]) {
        const auto a0 = ao_offsets[0] + ao[0];

        for(ao[1] = 0; ao[1] < n_aos1; ++ao[1]) {
            const auto a1 = ao_offsets[1] + ao[1];

            for(ao[2] = 0; ao[2] < n_aos2; ++ao[2]) {
                const auto a2 = ao_offsets[2] + ao[2];

                for(ao[3] = 0; ao[3] < n_aos3; ++ao[3]) {
                    const auto a3 = ao_offsets[3] + ao[3];

                    std::vector index{a0, a1, a2, a3};
                    t.set_elem(index, error);
                } // End loop over ao[3]
            } // End loop over ao[2]
        } // End loop over ao[1]
    } // End loop over ao[0]
}

} // namespace

using eri4_pt           = simde::ERI4;
using pair_estimator_pt = integrals::property_types::PrimitivePairEstimator;
using pt                = integrals::property_types::Uncertainty<eri4_pt>;

MODULE_CTOR(PrimitiveErrorModel) {
    satisfies_property_type<pt>();
    description(desc);
    // TODO citation for Chemist paper

    add_submodule<pair_estimator_pt>("Black Box Primitive Pair Estimator")
      .set_description("The module used to estimate the contributions of "
                       "primitive pairs to the overall integral values");
    add_submodule<pair_estimator_pt>("Primitive Pair Estimator")
      .set_description("The module used to estimate the contributions of "
                       "primitive pairs to the overall integral values");
}

MODULE_RUN(PrimitiveErrorModel) {
    const auto& [braket, tol] = pt::unwrap_inputs(inputs);

    const auto bra0 = braket.bra().first.ao_basis_set();
    const auto bra1 = braket.bra().second.ao_basis_set();
    const auto ket0 = braket.ket().first.ao_basis_set();
    const auto ket1 = braket.ket().second.ao_basis_set();

    // Get the Q_ab and Q_cd tensors (used to estimate error)
    auto& estimator     = submods.at("Primitive Pair Estimator");
    const auto& Q_ab_tw = estimator.run_as<pair_estimator_pt>(bra0, bra1);
    const auto& Q_cd_tw = estimator.run_as<pair_estimator_pt>(ket0, ket1);

    // Get the K_ab and K_cd tensors (used to determine which primitives are
    // neglected)
    auto& estimator2    = submods.at("Black Box Primitive Pair Estimator");
    const auto& K_ab_tw = estimator2.run_as<pair_estimator_pt>(bra0, bra1);
    const auto& K_cd_tw = estimator2.run_as<pair_estimator_pt>(ket0, ket1);

    // XXX: Workaround for needing the contiguous buffers to access elements
    using tensorwrapper::buffer::make_contiguous;
    const auto& Q_ab = make_contiguous(Q_ab_tw.buffer());
    const auto& Q_cd = make_contiguous(Q_cd_tw.buffer());
    const auto& K_ab = make_contiguous(K_ab_tw.buffer());
    const auto& K_cd = make_contiguous(K_cd_tw.buffer());

    // Get the number of AOs in each basis set
    const auto n_mu     = bra0.n_aos();
    const auto n_nu     = bra1.n_aos();
    const auto n_lambda = ket0.n_aos();
    const auto n_sigma  = ket1.n_aos();

    // Make the return buffer
    using float_type = double; // TODO: get from Q_ab or Q_cd
    tensorwrapper::shape::Smooth shape({n_mu, n_nu, n_lambda, n_sigma});
    std::vector<float_type> raw_buffer(shape.size(), 0.0);
    tensorwrapper::buffer::Contiguous buffer(std::move(raw_buffer), shape);

    // Make arrays to hold the indices and offsets
    using iter_type  = std::decay_t<decltype(n_mu)>; // Type of the loop index
    using shell_view = decltype(bra0.shell(0));
    std::array<iter_type, 4> n_shells{bra0.n_shells(), bra1.n_shells(),
                                      ket0.n_shells(), ket1.n_shells()};
    std::array<iter_type, 4> shell_i{0, 0, 0, 0};
    std::array<iter_type, 4> prim_offset{0, 0, 0, 0};
    std::array<iter_type, 4> ao_offsets{0, 0, 0, 0};

    // We'll collect the shells into this array
    std::array<shell_view, 4> shells;

    for(shell_i[0] = 0; shell_i[0] < n_shells[0]; ++shell_i[0]) {
        shells[0] = bra0.shell(shell_i[0]);

        prim_offset[1] = 0;
        ao_offsets[1]  = 0;
        for(shell_i[1] = 0; shell_i[1] < n_shells[1]; ++shell_i[1]) {
            shells[1] = bra1.shell(shell_i[1]);

            prim_offset[2] = 0;
            ao_offsets[2]  = 0;
            for(shell_i[2] = 0; shell_i[2] < n_shells[2]; ++shell_i[2]) {
                shells[2] = ket0.shell(shell_i[2]);

                prim_offset[3] = 0;
                ao_offsets[3]  = 0;
                for(shell_i[3] = 0; shell_i[3] < n_shells[3]; ++shell_i[3]) {
                    shells[3] = ket1.shell(shell_i[3]);

                    auto shell_error = shell_block_error(
                      shells, prim_offset, Q_ab, Q_cd, K_ab, K_cd, tol);
                    fill_ao_block(shells, ao_offsets, buffer, shell_error);

                    prim_offset[3] += shells[3].n_primitives();
                    ao_offsets[3] += shells[3].size();
                } // End loop over shell_i[3]

                prim_offset[2] += shells[2].n_primitives();
                ao_offsets[2] += shells[2].size();
            } // End loop over shell_i[2]

            prim_offset[1] += shells[1].n_primitives();
            ao_offsets[1] += shells[1].size();
        } // End loop over shell_i[1]

        prim_offset[0] += shells[0].n_primitives();
        ao_offsets[0] += shells[0].size();
    } // End loop over shell_i[0]

    simde::type::tensor error(shape, std::move(buffer));
    auto result = results();
    return pt::wrap_results(result, error);
}

} // namespace integrals::libint
#include "detail_/aos2shells.hpp"
#include "detail_/bases_helper.hpp"
#include "detail_/make_engine.hpp"
#include "detail_/make_shape.hpp"
#include "detail_/shells2ord.hpp"
#include "libint.hpp"
#include <simde/tensor_representation/ao_tensor_representation.hpp>

/// TODO: Unify implementations. Maybe with recursion?

namespace integrals {

using identity_op   = simde::type::el_identity;
using dipole_op     = simde::type::el_dipole;
using quadrupole_op = simde::type::el_quadrupole;
using octupole_op   = simde::type::el_octupole;

using overlap_pt    = simde::AOTensorRepresentation<2, identity_op>;
using dipole_pt     = simde::AOTensorRepresentation<2, dipole_op>;
using quadrupole_pt = simde::AOTensorRepresentation<2, quadrupole_op>;
using octupole_pt   = simde::AOTensorRepresentation<2, octupole_op>;

/// Grab the various detail_ functions
using namespace detail_;

template<std::size_t L, typename OperatorType>
TEMPLATED_MODULE_CTOR(LibintMultipole, L, OperatorType) {
    description("Computes an in-core integral with libint");

    /// This should satisfy overlap, but we can't reduce the dimensionality
    /// of the tensor at the moment.
    // satisfies_property_type<overlap_pt>();
    // identity_op I;
    // change_input(I.as_string()).change(std::move(I));

    satisfies_property_type<dipole_pt>();
    dipole_op r;
    change_input(r.as_string()).change(std::move(r));

    if constexpr(L > 0) {
        satisfies_property_type<quadrupole_pt>();
        quadrupole_op r2;
        change_input(r2.as_string()).change(std::move(r2));
    }

    if constexpr(L > 1) {
        satisfies_property_type<octupole_pt>();
        octupole_op r3;
        change_input(r3.as_string()).change(std::move(r3));
    }

    add_input<double>("Threshold")
      .set_default(1.0E-16)
      .set_description(
        "The target precision with which the integrals will be computed");
}

template<std::size_t L, typename OperatorType>
TEMPLATED_MODULE_RUN(LibintMultipole, L, OperatorType) {
    /// Typedefs
    using size_vector_t = std::vector<std::size_t>;
    using tensor_t      = simde::type::tensor;
    using field_t       = typename tensor_t::field_type;

    /// Grab input information
    auto bases  = unpack_bases<2>(inputs);
    auto op_str = OperatorType().as_string();
    auto op     = inputs.at(op_str).template value<const OperatorType&>();
    auto thresh = inputs.at("Threshold").value<double>();

    /// Lambda to calculate values
    auto l = [&](const auto& lo, const auto& up, auto* data) {
        /// Convert index values from AOs to shells
        /// Leading index is for multipole components
        constexpr std::size_t N = 2;
        size_vector_t lo_shells, up_shells;
        for(auto i = 0; i < N; ++i) {
            auto shells_in_tile = aos2shells(bases[i], lo[i + 1], up[i + 1]);
            lo_shells.push_back(shells_in_tile.front());
            up_shells.push_back(shells_in_tile.back());
        }
        // for(auto& i : lo_shells) std::cout << i << " ";
        // std::cout << std::endl;
        // for(auto& i : up_shells) std::cout << i << " ";
        // std::cout << std::endl;

        /// Calculate the number of values per leading index
        auto leading_step = 0;
        for(auto i = lo_shells[0]; i <= up_shells[0]; ++i) {
            for(auto j = lo_shells[1]; j <= up_shells[1]; ++j) {
                leading_step += bases[0][i].size() * bases[1][j].size();
            }
        }
        /// std::cout << leading_step << std::endl;

        /// Make the libint engine to calculate integrals
        auto engine     = make_engine(bases, op, thresh);
        const auto& buf = engine.results();

        /// Loop through shell combinations
        size_vector_t curr_shells = lo_shells;
        while(curr_shells[0] <= up_shells[0]) {
            /// Determine which values will be computed this time
            auto ord_pos = shells2ord(bases, curr_shells, lo_shells, up_shells);

            /// Compute values
            engine.compute(bases[0][curr_shells[0]], bases[1][curr_shells[1]]);

            /// Copy libint values into tile data;
            for(auto i = lo[0]; i < up[0]; ++i) {
                auto depth = i * leading_step;
                for(auto j = 0; j < ord_pos.size(); ++j) {
                    data[ord_pos[j] + depth] = buf[i][j];
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

    /// Count up necessary components for multipole
    std::size_t leading_extent = 4;
    if constexpr(L == 1) { leading_extent = 10; }
    if constexpr(L == 2) { leading_extent = 20; }

    /// Make complete tensor and slice out return values
    tensor_t I(l, make_shape(bases, leading_extent),
               tensorwrapper::tensor::default_allocator<field_t>());

    auto rv = results();
    auto D  = I.slice({1, 0, 0}, {4, 7, 7});
    rv      = dipole_pt::wrap_results(rv, D);
    if constexpr(L > 0) {
        auto Q = I.slice({4, 0, 0}, {10, 7, 7});
        rv     = quadrupole_pt::wrap_results(rv, Q);
    }
    if constexpr(L > 1) {
        auto O = I.slice({10, 0, 0}, {20, 7, 7});
        rv     = octupole_pt::wrap_results(rv, O);
    }
    return rv;
}

template class LibintMultipole<0, simde::type::el_dipole>;
template class LibintMultipole<1, simde::type::el_quadrupole>;
template class LibintMultipole<2, simde::type::el_octupole>;

} // namespace integrals

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

MODULE_CTOR(LibintDipole) {
    description("Computes in-core dipole integrals with libint");

    /// This should satisfy overlap, but we can't reduce the dimensionality
    /// of the tensor at the moment.
    // satisfies_property_type<overlap_pt>();
    // identity_op I;
    // change_input(I.as_string()).change(std::move(I));

    satisfies_property_type<dipole_pt>();
    dipole_op r;
    change_input(r.as_string()).change(std::move(r));

    add_input<double>("Threshold")
      .set_default(1.0E-16)
      .set_description(
        "The target precision with which the integrals will be computed");
}

MODULE_RUN(LibintDipole) {
    /// Typedefs
    using size_vector_t = std::vector<std::size_t>;
    using tensor_t      = simde::type::tensor;
    using field_t       = typename tensor_t::field_type;

    /// Grab input information
    auto bases  = unpack_bases<2>(inputs);
    auto op_str = dipole_op().as_string();
    auto op     = inputs.at(op_str).template value<const dipole_op&>();
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
        ///std::cout << leading_step << std::endl;

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
    tensor_t I(l, make_shape(bases, std::size_t{4}),
               tensorwrapper::tensor::default_allocator<field_t>());
    auto D = I.slice({1, 0, 0}, {4, 7, 7});

    auto rv = results();
    rv      = dipole_pt::wrap_results(rv, D);
    return rv;
}

MODULE_CTOR(LibintQuadrupole) {
    description("Computes an in-core integral with libint");

    /// This should satisfy overlap, but we can't reduce the dimensionality
    /// of the tensor at the moment.
    // satisfies_property_type<overlap_pt>();
    // identity_op I;
    // change_input(I.as_string()).change(std::move(I));

    satisfies_property_type<dipole_pt>();
    dipole_op r;
    change_input(r.as_string()).change(std::move(r));

    satisfies_property_type<quadrupole_pt>();
    quadrupole_op r2;
    change_input(r2.as_string()).change(std::move(r2));

    add_input<double>("Threshold")
      .set_default(1.0E-16)
      .set_description(
        "The target precision with which the integrals will be computed");
}

MODULE_RUN(LibintQuadrupole) {
    /// Typedefs
    using size_vector_t = std::vector<std::size_t>;
    using tensor_t      = simde::type::tensor;
    using field_t       = typename tensor_t::field_type;

    /// Grab input information
    auto bases  = unpack_bases<2>(inputs);
    auto op_str = quadrupole_op().as_string();
    auto op     = inputs.at(op_str).template value<const quadrupole_op&>();
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
        ///std::cout << leading_step << std::endl;

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
    tensor_t I(l, make_shape(bases, std::size_t{10}),
               tensorwrapper::tensor::default_allocator<field_t>());
    auto D = I.slice({1, 0, 0}, {4, 7, 7});
    auto Q = I.slice({4, 0, 0}, {10, 7, 7});

    auto rv = results();
    rv      = dipole_pt::wrap_results(rv, D);
    rv      = quadrupole_pt::wrap_results(rv, Q);
    return rv;
}

MODULE_CTOR(LibintOctupole) {
    description("Computes an in-core integral with libint");

    /// This should satisfy overlap, but we can't reduce the dimensionality
    /// of the tensor at the moment.
    // satisfies_property_type<overlap_pt>();
    // identity_op I;
    // change_input(I.as_string()).change(std::move(I));

    satisfies_property_type<dipole_pt>();
    dipole_op r;
    change_input(r.as_string()).change(std::move(r));

    satisfies_property_type<quadrupole_pt>();
    quadrupole_op r2;
    change_input(r2.as_string()).change(std::move(r2));

    satisfies_property_type<octupole_pt>();
    octupole_op r3;
    change_input(r3.as_string()).change(std::move(r3));

    add_input<double>("Threshold")
      .set_default(1.0E-16)
      .set_description(
        "The target precision with which the integrals will be computed");
}

MODULE_RUN(LibintOctupole) {
    /// Typedefs
    using size_vector_t = std::vector<std::size_t>;
    using tensor_t      = simde::type::tensor;
    using field_t       = typename tensor_t::field_type;

    /// Grab input information
    auto bases  = unpack_bases<2>(inputs);
    auto op_str = octupole_op().as_string();
    auto op     = inputs.at(op_str).template value<const octupole_op&>();
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
        ///std::cout << leading_step << std::endl;

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
    tensor_t I(l, make_shape(bases, std::size_t{20}),
               tensorwrapper::tensor::default_allocator<field_t>());
    auto D = I.slice({1, 0, 0}, {4, 7, 7});
    auto Q = I.slice({4, 0, 0}, {10, 7, 7});
    auto O = I.slice({10, 0, 0}, {20, 7, 7});

    auto rv = results();
    rv      = dipole_pt::wrap_results(rv, D);
    rv      = quadrupole_pt::wrap_results(rv, Q);
    rv      = octupole_pt::wrap_results(rv, O);
    return rv;
}

} // namespace integrals

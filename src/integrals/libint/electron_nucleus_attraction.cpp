#include "detail_/aos2shells.hpp"
#include "detail_/bases_helper.hpp"
#include "detail_/make_engine.hpp"
#include "detail_/make_shape.hpp"
#include "detail_/shells2ord.hpp"
#include "libint.hpp"
#include <simde/tensor_representation/ao_tensor_representation.hpp>

namespace integrals {

/// Grab the various detail_ functions
using namespace detail_;

MODULE_CTOR(LibintDOI) {
    description("Computes DOI integrals with Libint");
    using my_pt = simde::AOTensorRepresentation<2, simde::type::el_el_delta>;

    satisfies_property_type<my_pt>();

    add_input<double>("Threshold")
      .set_default(1.0E-16)
      .set_description(
        "The target precision with which the integrals will be computed");
}

MODULE_RUN(LibintDOI) {
    using op_t          = simde::type::el_el_delta;
    using my_pt         = simde::AOTensorRepresentation<2, op_t>;
    using size_vector_t = std::vector<std::size_t>;
    using tensor_t      = simde::type::tensor;
    using field_t       = typename tensor_t::field_type;

    auto init_bases = unpack_bases<2>(inputs);
    auto op_str     = op_t().as_string();
    auto op         = inputs.at(op_str).template value<const op_t&>();
    auto thresh     = inputs.at("Threshold").value<double>();

    /// Have to double up the basis sets
    std::vector<libint2::BasisSet> bases;
    for(auto& set : init_bases) {
        for(auto i = 0; i < 2; ++i) bases.push_back(set);
    }

    /// Lambda to calculate values
    auto l = [&](const auto& lo, const auto& up, auto* data) {
        /// Convert index values from AOs to shells
        constexpr std::size_t N = 4;
        size_vector_t lo_shells, up_shells;
        for(auto i = 0; i < N; ++i) {
            auto shells_in_tile = aos2shells(bases[i], lo[i], up[i]);
            lo_shells.push_back(shells_in_tile.front());
            up_shells.push_back(shells_in_tile.back());
        }

        /// Make the libint engine to calculate integrals
        auto engine     = make_engine(bases, op, thresh);
        const auto& buf = engine.results();

        /// Loop through shell combinations
        size_vector_t curr_shells = lo_shells;
        while(curr_shells[0] <= up_shells[0]) {
            /// Determine which values will be computed this time
            auto ord_pos = shells2ord(bases, curr_shells);

            /// Compute values
            run_engine_(engine, bases, curr_shells,
                        std::make_index_sequence<N>());
            auto vals = buf[0];

            /// Copy libint values into tile data;
            for(auto i = 0; i < ord_pos.size(); ++i) {
                data[ord_pos[i]] = vals[i];
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
    tensor_t I(l, make_shape(bases),
               tensorwrapper::tensor::default_allocator<field_t>());

    /// Finish
    auto rv = results();
    return my_pt::wrap_results(rv, I);
}

} // namespace integrals

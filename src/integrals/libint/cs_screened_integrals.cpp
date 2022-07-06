#include "cs_screened_integrals.hpp"
#include "detail_/aos2shells.hpp"
#include "detail_/bases_helper.hpp"
#include "detail_/make_engine.hpp"
#include "detail_/make_shape.hpp"
#include "detail_/shells2ord.hpp"
#include <simde/cauchy_schwarz_approximation.hpp>
#include <simde/tensor_representation/ao_tensor_representation.hpp>

namespace integrals {

/// Grab the various detail_ functions
using namespace detail_;

template<std::size_t N, typename OperatorType>
TEMPLATED_MODULE_CTOR(CSLibint, N, OperatorType) {
    description("Computes an in-core integral with libint");
    using my_pt = simde::AOTensorRepresentation<N, OperatorType>;

    satisfies_property_type<my_pt>();

    add_input<double>("Threshold")
      .set_default(1.0E-16)
      .set_description(
        "The target precision with which the integrals will be computed");

    add_input<double>("Screening Threshold")
      .set_default(0.0)
      .set_description("Cauchy-Schwarz Screening Threshold");

    using cs_approx_pt = simde::ShellNorms<OperatorType>;
    add_submodule<cs_approx_pt>("Shell Norms")
      .set_description(
        "Computes the Cauchy-Schwarz Matrix for a pair of basis sets");
}

template<std::size_t N, typename OperatorType>
TEMPLATED_MODULE_RUN(CSLibint, N, OperatorType) {
    /// Typedefs
    using my_pt         = simde::AOTensorRepresentation<N, OperatorType>;
    using size_vector_t = std::vector<std::size_t>;
    using tensor_t      = simde::type::tensor;
    using field_t       = typename tensor_t::field_type;

    /// Grab input information
    auto bases     = unpack_bases<N>(inputs);
    auto op_str    = OperatorType().as_string();
    auto op        = inputs.at(op_str).template value<const OperatorType&>();
    auto thresh    = inputs.at("Threshold").value<double>();
    auto cs_thresh = inputs.at("Screening Threshold").value<double>();

    /// Calculate Shell Norms for screening
    auto& cs_screen = submods.at("Shell Norms");

    /// Lambda to calculate values
    auto l = [&](const auto& lo, const auto& up, auto* data) {
        /// Convert index values from AOs to shells
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
            /// Check if current shells screen out

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

    /// Geminal exponent handling
    constexpr auto is_stg =
      std::is_same_v<OperatorType, simde::type::el_el_stg>;
    constexpr auto is_yukawa =
      std::is_same_v<OperatorType, simde::type::el_el_yukawa>;
    if constexpr(is_stg || is_yukawa) {
        auto I_ann = I(I.make_annotation());
        I_ann      = op.template at<0>().coefficient * I_ann;
    }

    /// Finish
    auto rv = results();
    return my_pt::wrap_results(rv, I);
}

template class CSLibint<3, simde::type::el_el_coulomb>;
template class CSLibint<4, simde::type::el_el_coulomb>;
template class CSLibint<3, simde::type::el_el_stg>;
template class CSLibint<4, simde::type::el_el_stg>;
template class CSLibint<3, simde::type::el_el_yukawa>;
template class CSLibint<4, simde::type::el_el_yukawa>;

} // namespace integrals

#include "detail_/aos2shells.hpp"
#include "detail_/bases_helper.hpp"
#include "detail_/make_engine.hpp"
#include "detail_/make_shape.hpp"
#include "detail_/shells2ord.hpp"
#include "libint.hpp"
#include <simde/tensor_representation/ao_tensor_representation.hpp>

namespace integrals {

using size_vector_t  = std::vector<std::size_t>;
using bases_vector_t = std::vector<libint2::BasisSet>;
using tensor_t       = simde::type::tensor;
using field_t        = typename tensor_t::field_type;

/** @brief Wrap the call of LibInt2 engine so it can take a variable number
 * of shell inputs.
 *
 * @tparam Is A variadic parameter pack of integers from [0,NBases) to
 * expand.
 * @param engine The LibInt2 engine that computes integrals
 * @param bases The bases sets that hold the shells
 * @param shells The index of the requested shell block
 */
template<std::size_t... Is>
void run_engine_(libint2::Engine& engine, bases_vector_t& bases,
                 size_vector_t& shells, std::index_sequence<Is...>) {
    engine.compute(bases[Is][shells[Is]]...);
}

template<std::size_t N, typename OperatorType>
TEMPLATED_MODULE_CTOR(Libint, N, OperatorType) {
    description("Computes integrals with Libint");
    using my_pt = simde::AOTensorRepresentation<N, OperatorType>;

    satisfies_property_type<my_pt>();

    add_input<double>("Threshold")
      .set_default(1.0E-16)
      .set_description(
        "The target precision with which the integrals will be computed");
}

template<std::size_t N, typename OperatorType>
TEMPLATED_MODULE_RUN(Libint, N, OperatorType) {
    using my_pt = simde::AOTensorRepresentation<N, OperatorType>;

    auto bases  = detail_::unpack_bases<N>(inputs);
    auto op_str = OperatorType().as_string();
    auto op     = inputs.at(op_str).template value<const OperatorType&>();
    auto thresh = inputs.at("Threshold").value<double>();

    auto l = [&](const auto& lo, const auto& up, auto* data) {
        /// Convert index values from AOs to shells
        size_vector_t lo_shells, up_shells;
        for(auto i = 0; i < N; ++i) {
            auto shells_in_tile = detail_::aos2shells(bases[i], lo[i], up[i]);
            lo_shells.push_back(shells_in_tile.front());
            up_shells.push_back(shells_in_tile.back());
        }

        /// Make the libint engine to calculate integrals
        auto engine     = detail_::make_engine(bases, op, thresh);
        const auto& buf = engine.results();

        /// Loop through shell combinations
        size_vector_t curr_shells = lo_shells;
        while(curr_shells[0] <= up_shells[0]) {
            /// Determine which values will be computed this time
            auto ord_pos = detail_::shells2ord(bases, curr_shells);

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
    tensor_t I(l, detail_::make_shape(bases),
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

    auto rv = results();
    return my_pt::wrap_results(rv, I);
}

template class Libint<2, simde::type::el_el_coulomb>;
template class Libint<3, simde::type::el_el_coulomb>;
template class Libint<4, simde::type::el_el_coulomb>;
template class Libint<2, simde::type::el_kinetic>;
template class Libint<2, simde::type::el_nuc_coulomb>;
template class Libint<2, simde::type::el_identity>;
template class Libint<2, simde::type::el_el_stg>;
template class Libint<3, simde::type::el_el_stg>;
template class Libint<4, simde::type::el_el_stg>;
template class Libint<2, simde::type::el_el_yukawa>;
template class Libint<3, simde::type::el_el_yukawa>;
template class Libint<4, simde::type::el_el_yukawa>;
template class Libint<2, simde::type::el_el_f12_commutator>;
template class Libint<3, simde::type::el_el_f12_commutator>;
template class Libint<4, simde::type::el_el_f12_commutator>;

} // namespace integrals

#include "detail_/bases_helper.hpp"
#include "libint.hpp"
#include <simde/tensor_representation/ao_tensor_representation.hpp>

namespace integrals {

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

    /// Fill in here

    simde::type::tensor I;

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
template class Libint<4, simde::type::el_el_f12_commutator>;
} // namespace integrals

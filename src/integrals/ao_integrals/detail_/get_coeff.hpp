#pragma once
#include <simde/types.hpp>
#include <type_traits>

namespace integrals::ao_integrals::detail_ {

template<typename OperatorType>
auto get_coefficient(const OperatorType& op) {
    /// Geminal exponent handling
    constexpr auto is_stg =
      std::is_same_v<OperatorType, simde::type::el_el_stg>;
    constexpr auto is_yukawa =
      std::is_same_v<OperatorType, simde::type::el_el_yukawa>;

    double coeff = 1.0;
    if constexpr(is_stg || is_yukawa) {
        coeff = op.template at<0>().coefficient;
    }
    return coeff;
}

} // namespace integrals::ao_integrals

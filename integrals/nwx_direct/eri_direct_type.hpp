#pragma once
#include <property_types/ao_integrals/electron_repulsion.hpp>
#include <sde/property_type/property_type.hpp>

namespace property_types {

template<typename BaseType>
DECLARE_DERIVED_TEMPLATED_PROPERTY_TYPE(Direct, BaseType, BaseType);

template<typename BaseType>
TEMPLATED_PROPERTY_TYPE_INPUTS(Direct, BaseType) {
    return sde::declare_input();
}

template<typename BaseType>
TEMPLATED_PROPERTY_TYPE_RESULTS(Direct, BaseType) {
    return sde::declare_result();
}

extern template class Direct<ao_integrals::ERI2C<float>>;
extern template class Direct<ao_integrals::ERI2C<double>>;
extern template class Direct<ao_integrals::ERI3C<float>>;
extern template class Direct<ao_integrals::ERI3C<double>>;
extern template class Direct<ao_integrals::ERI4C<float>>;
extern template class Direct<ao_integrals::ERI4C<double>>;

} // namespace property_types

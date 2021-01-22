#pragma once
#include "integrals/types.hpp"
#include <property_types/ao_integrals/electron_repulsion.hpp>
#include <sde/module_base.hpp>

namespace integrals {

template<typename BaseType>
DECLARE_MODULE(LibintDirect);

template<typename ElementType>
using DirectERI3C =
  LibintDirect<property_types::ao_integrals::ERI3C<ElementType>>;

template<typename ElementType>
using DirectERI4C =
  LibintDirect<property_types::ao_integrals::ERI4C<ElementType>>;

extern template class LibintDirect<property_types::ao_integrals::ERI3C<double>>;
extern template class LibintDirect<property_types::ao_integrals::ERI4C<double>>;

using ERI3Direct = DirectERI3C<double>;
using ERI4Direct = DirectERI4C<double>;
} // namespace integrals
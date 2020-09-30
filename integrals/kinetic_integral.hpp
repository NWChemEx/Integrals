#pragma once
#include "integrals/types.hpp"
#include <sde/module_base.hpp>

namespace integrals {

template<typename element_type = double>
DECLARE_MODULE(KineticInt);

extern template class KineticInt<double>;

using Kinetic = KineticInt<double>;
} // namespace integrals
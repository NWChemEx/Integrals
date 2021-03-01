#pragma once
#include "integrals/types.hpp"
#include <sde/module_base.hpp>

namespace integrals {

template<typename element_type = double>
DECLARE_MODULE(NuclearInt);

extern template class NuclearInt<double>;

using Nuclear = NuclearInt<double>;
} // namespace integrals
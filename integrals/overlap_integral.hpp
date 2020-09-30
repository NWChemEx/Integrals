#pragma once
#include "integrals/types.hpp"
#include <sde/module_base.hpp>

namespace integrals {

template<typename element_type = double>
DECLARE_MODULE(OverlapInt);

extern template class OverlapInt<double>;

using Overlap = OverlapInt<double>;
} // namespace integrals
#pragma once
#include "integrals/types.hpp"
#include <sde/module_base.hpp>

namespace integrals {

template<typename element_type = double>
DECLARE_MODULE(DOInt);

extern template class DOInt<double>;

using DOI = DOInt<double>;
} // namespace integrals
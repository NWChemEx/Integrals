#pragma once
#include "integrals/types.hpp"
#include <sde/module_base.hpp>

namespace integrals {

template<typename element_type = double>
DECLARE_MODULE(STG3CIntDirect);

template<typename element_type = double>
DECLARE_MODULE(STG4CIntDirect);

extern template class STG3CIntDirect<double>;
extern template class STG4CIntDirect<double>;

using STG3Direct = STG3CIntDirect<double>;
using STG4Direct = STG4CIntDirect<double>;
} // namespace integrals
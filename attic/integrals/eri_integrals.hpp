#pragma once
#include <sde/module_base.hpp>

namespace integrals {

template<typename element_type = double>
DECLARE_MODULE(ERI2CInt);

template<typename element_type = double>
DECLARE_MODULE(ERI3CInt);

template<typename element_type = double>
DECLARE_MODULE(ERI4CInt);

extern template class ERI2CInt<double>;
extern template class ERI3CInt<double>;
extern template class ERI4CInt<double>;

using ERI2 = ERI2CInt<double>;
using ERI3 = ERI3CInt<double>;
using ERI4 = ERI4CInt<double>;
} // namespace integrals
#pragma once
#include <sde/module_base.hpp>

namespace integrals {

template<typename element_type = double>
DECLARE_MODULE(Yukawa3CIntDirect);

template<typename element_type = double>
DECLARE_MODULE(Yukawa4CIntDirect);

extern template class Yukawa3CIntDirect<double>;
extern template class Yukawa4CIntDirect<double>;

using Yukawa3Direct = Yukawa3CIntDirect<double>;
using Yukawa4Direct = Yukawa4CIntDirect<double>;
} // namespace integrals
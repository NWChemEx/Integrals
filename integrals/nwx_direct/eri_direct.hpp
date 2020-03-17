#pragma once
#include <sde/module_base.hpp>
#include "integrals/types.hpp"

namespace integrals {

    template<typename element_type = double>
    DECLARE_MODULE(ERI3CIntDirect);

    template<typename element_type = double>
    DECLARE_MODULE(ERI4CIntDirect);

    extern template class ERI3CIntDirect<double>;
    extern template class ERI4CIntDirect<double>;

    using ERI3Direct = ERI3CIntDirect<double>;
    using ERI4Direct = ERI4CIntDirect<double>;
} // namespace integrals
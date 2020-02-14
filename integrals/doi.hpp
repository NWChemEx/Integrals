#pragma once
#include <sde/module_base.hpp>
#include "integrals/types.hpp"

namespace integrals {

    template<typename element_type = double>
    DECLARE_MODULE(DOInt);

    extern template class DOInt<double>;

    using DOI = DOInt<double>;
} // namespace integrals
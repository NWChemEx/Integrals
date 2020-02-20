#pragma once
#include <sde/module_base.hpp>
#include "integrals/types.hpp"

namespace integrals {

    template<typename element_type = double>
    DECLARE_MODULE(KineticInt);

    extern template class KineticInt<double>;

    using Kinetic = KineticInt<double>;
} // namespace integrals
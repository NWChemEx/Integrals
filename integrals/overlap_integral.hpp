#pragma once
#include <sde/module_base.hpp>
#include <property_types/ao_integrals/overlap.hpp>
#include "integrals/types.hpp"

namespace integrals {

    template<typename element_type = double>
    DECLARE_MODULE(OverlapInt);

    extern template class OverlapInt<double>;

    using Overlap = OverlapInt<double>;
} // namespace integrals
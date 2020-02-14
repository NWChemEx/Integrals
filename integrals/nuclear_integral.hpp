#pragma once
#include <sde/module_base.hpp>
#include <property_types/ao_integrals/nuclear.hpp>
#include "integrals/types.hpp"

namespace integrals {

    template<typename element_type = double>
    DECLARE_MODULE(NuclearInt);

    extern template class NuclearInt<double>;

    using Nuclear = NuclearInt<double>;
} // namespace integrals
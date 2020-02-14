#pragma once
#include <sde/module_base.hpp>
#include <property_types/ao_integrals/stg.hpp>
#include "integrals/types.hpp"

namespace integrals {

    template<typename element_type = double>
    DECLARE_MODULE(STG2CInt);

    template<typename element_type = double>
    DECLARE_MODULE(STG3CInt);

    template<typename element_type = double>
    DECLARE_MODULE(STG4CInt);

    extern template class STG2CInt<double>;
    extern template class STG3CInt<double>;
    extern template class STG4CInt<double>;

    using STG2 = STG2CInt<double>;
    using STG3 = STG3CInt<double>;
    using STG4 = STG4CInt<double>;
} // namespace integrals
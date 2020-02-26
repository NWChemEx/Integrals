#pragma once
#include <sde/module_base.hpp>
#include "integrals/types.hpp"

namespace integrals {

    template<typename element_type = double>
    DECLARE_MODULE(Yukawa2CInt);

    template<typename element_type = double>
    DECLARE_MODULE(Yukawa3CInt);

    template<typename element_type = double>
    DECLARE_MODULE(Yukawa4CInt);

    extern template class Yukawa2CInt<double>;
    extern template class Yukawa3CInt<double>;
    extern template class Yukawa4CInt<double>;

    using Yukawa2 = Yukawa2CInt<double>;
    using Yukawa3 = Yukawa3CInt<double>;
    using Yukawa4 = Yukawa4CInt<double>;
} // namespace integrals
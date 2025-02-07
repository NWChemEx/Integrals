#pragma once

#ifdef ENABLE_SIGMA
#include <sigma/sigma.hpp>
#endif

namespace integrals::type {
#ifdef ENABLE_SIGMA
using uncertain_float  = sigma::UFloat;
using uncertain_double = sigma::UDouble;
#else
using uncertain_float  = void;
using uncertain_double = void;
#endif

} // namespace integrals::type
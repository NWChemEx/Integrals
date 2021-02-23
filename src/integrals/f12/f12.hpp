#pragma once
#include "integrals/property_types.hpp"
#include <sde/module_base.hpp>

#define DECLARE_MULTICENTER_STG(TYPEDEF, PROP_TYPE)                    \
    template<typename PropType>                                        \
    DECLARE_MODULE(STG##PROP_TYPE);                                    \
    template<typename T>                                               \
    using stg_##TYPEDEF##2c = STG##PROP_TYPE<pt::TYPEDEF##2c < T> > ;  \
    template<typename T>                                               \
    using stg_##TYPEDEF##3c = STG##PROP_TYPE<pt::TYPEDEF##3c < T> > ;  \
    template<typename T>                                               \
    using stg_##TYPEDEF##4c = STG##PROP_TYPE<pt::TYPEDEF##4c < T> > ;  \
    extern template class STG##PROP_TYPE<pt::TYPEDEF##2c < double> > ; \
    extern template class STG##PROP_TYPE<pt::TYPEDEF##3c < double> > ; \
    extern template class STG##PROP_TYPE<pt::TYPEDEF##4c < double> >

/** @brief Namespace for derived integral types needed in explicitly correlated
 *         theories.
 */
namespace integrals::f12 {

DECLARE_MULTICENTER_STG(correlation_factor_, CorrelationFactor);
DECLARE_MULTICENTER_STG(gr, GR);

} // namespace integrals::f12

#undef DECLARE_MULTICENTER_STG
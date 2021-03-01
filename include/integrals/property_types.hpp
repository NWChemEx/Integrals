/** @file property_types.hpp
 *
 *  This file establishes typedefs for the property types used throughout the
 *  Integrals repos. For example we create the alias `pt::overlap<T>` for the
 *  type of the overlap matrix in the AO basis set. Files in this repo should
 *  use the typedefs here to avoid coupling directly to the property types.
 */

#pragma once
#include <property_types/ao_integrals/ao_integrals.hpp>
#include <property_types/transformed.hpp>

/** @brief Code factorization for declaring a typedef of a property type that is
 *         defined for 2, 3, and 4 centers.
 *
 *  This macro replaces the boilerplate associated with declaring typdefs for
 *  property types that are defined for a number of centers.
 *
 *  For example:
 *  ```.cpp
 *  MULTICENTER_TYPEDEF(eri, ERI);
 *  ```
 *  defines three templated typedefs: `eri2c`, `eri3c`, and `eri4c` which are
 *  typedefs of `property_types::ao_integrals`: `ERI2c`, `ERI3c`, and `ERI4c`
 *  classes.
 *
 *  @param[in] TYPEDEF The prefix for the typedef to define. The resulting
 *                     typedefs will be `TYPEDEF2c`, `TYPEDEF3c`, and
 *                     `TYPEDEF4c`.
 *  @param[in] PROP_TYPE The prefix for class defining the property type. The
 *                       class is assumed to live in the
 *                       `property_types::ao_integrals` namespace and have a
 *                       2, 3, and 4 center version defined.
 */
#define MULTICENTER_TYPEDEF(TYPEDEF, PROP_TYPE)                     \
    template<typename ElementType>                                  \
    using TYPEDEF##2c =                                             \
      property_types::ao_integrals::PROP_TYPE##2C < ElementType > ; \
    template<typename ElementType>                                  \
    using TYPEDEF##3c =                                             \
      property_types::ao_integrals::PROP_TYPE##3C < ElementType > ; \
    template<typename ElementType>                                  \
    using TYPEDEF##4c =                                             \
      property_types::ao_integrals::PROP_TYPE##4C < ElementType >

namespace integrals::pt {

/// Property types for the F12 quantity f12
MULTICENTER_TYPEDEF(correlation_factor_, f12::CorrelationFactor);

/// Property types for the F12 quantity f12*f12
MULTICENTER_TYPEDEF(correlation_factor_squared_, f12::CorrelationFactorSquared);

/// Property types for the F12 quantity [f12, [T, f12]]
MULTICENTER_TYPEDEF(dfdr_squared_, f12::dfdrSquared);

/// Property type for the differential overlap integral
template<typename ElementType>
using doi = property_types::ao_integrals::DOI<ElementType>;

/// Property type for the electronic dipole moment
template<typename ElementType>
using edipole = property_types::ao_integrals::EDipole<ElementType>;

/// Property type for the electronic quadrupole moment
template<typename ElementType>
using equadrupole = property_types::ao_integrals::EQuadrupole<ElementType>;

/// Property type for the electronic octopole moment
template<typename ElementType>
using eoctopole = property_types::ao_integrals::EOctopole<ElementType>;

/// Property types for electron repulsion integrals
MULTICENTER_TYPEDEF(eri, ERI);

/// Property types for the F12 quantity f12/r12
MULTICENTER_TYPEDEF(gr, f12::GR);

/// Property type for the electronic kinetic energy
template<typename ElementType>
using kinetic = property_types::ao_integrals::Kinetic<ElementType>;

/// Property type for the electron-nuclear attraction
template<typename ElementType>
using nuclear = property_types::ao_integrals::Nuclear<ElementType>;

template<typename ElementType>
using overlap = property_types::ao_integrals::Overlap<ElementType>;

MULTICENTER_TYPEDEF(stg, STG);

template<typename BaseType>
using transformed = property_types::Transformed<BaseType>;

MULTICENTER_TYPEDEF(yukawa, Yukawa);

} // namespace integrals::pt

#undef MULTICENTER_TYPEDEF
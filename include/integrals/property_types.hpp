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

namespace integrals::pt {

template<typename ElementType>
using doi = property_types::ao_integrals::DOI<ElementType>;

template<typename ElementType>
using edipole = property_types::ao_integrals::EDipole<ElementType>;

template<typename ElementType>
using equadrupole = property_types::ao_integrals::EQuadrupole<ElementType>;

template<typename ElementType>
using eoctopole = property_types::ao_integrals::EOctopole<ElementType>;

template<typename ElementType>
using eri2c = property_types::ao_integrals::ERI2C<ElementType>;

template<typename ElementType>
using eri3c = property_types::ao_integrals::ERI3C<ElementType>;

template<typename ElementType>
using eri4c = property_types::ao_integrals::ERI4C<ElementType>;

template<typename ElementType>
using kinetic = property_types::ao_integrals::Kinetic<ElementType>;

template<typename ElementType>
using nuclear = property_types::ao_integrals::Nuclear<ElementType>;

template<typename ElementType>
using overlap = property_types::ao_integrals::Overlap<ElementType>;

template<typename ElementType>
using stg2c = property_types::ao_integrals::STG2C<ElementType>;

template<typename ElementType>
using stg3c = property_types::ao_integrals::STG3C<ElementType>;

template<typename ElementType>
using stg4c = property_types::ao_integrals::STG4C<ElementType>;

template<typename BaseType>
using transformed = property_types::Transformed<BaseType>;

template<typename ElementType>
using yukawa2c = property_types::ao_integrals::Yukawa2C<ElementType>;

template<typename ElementType>
using yukawa3c = property_types::ao_integrals::Yukawa3C<ElementType>;

template<typename ElementType>
using yukawa4c = property_types::ao_integrals::Yukawa4C<ElementType>;

} // namespace integrals::pt
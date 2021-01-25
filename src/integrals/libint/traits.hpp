#pragma once
#include "integrals/property_types.hpp"
#include <type_traits>

// TODO: The non-libint traits should probably be in property types to help
// other people who want to do meta-template programming with property types

namespace integrals {
namespace detail_ {

template<typename T>
struct LibintOp;

template<typename T>
struct LibintOp<pt::doi<T>> {
    static constexpr auto value = libint2::Operator::delta;
};

template<typename T>
struct LibintOp<property_types::ao_integrals::ERI<T>> {
    static constexpr auto value = libint2::Operator::coulomb;
};

template<typename T>
struct LibintOp<property_types::ao_integrals::Kinetic<T>> {
    static constexpr auto value = libint2::Operator::kinetic;
};

template<typename T>
struct LibintOp<property_types::ao_integrals::Overlap<T>> {
    static constexpr auto value = libint2::Operator::overlap;
};

template<typename T>
struct LibintOp<property_types::ao_integrals::Nuclear<T>> {
    static constexpr auto value = libint2::Operator::nuclear;
};

template<typename T>
struct LibintOp<property_types::ao_integrals::STG<T>> {
    static constexpr auto value = libint2::Operator::stg;
};

template<typename T>
struct LibintOp<property_types::ao_integrals::Yukawa<T>> {
    static constexpr auto value = libint2::Operator::yukawa;
};

template<typename T>
struct IsDOI : std::false_type {};

template<typename T>
struct IsDOI<pt::doi<T>> : std::true_type {};

template<typename T>
struct IsNuclear : std::false_type {};

template<typename T>
struct IsNuclear<pt::nuclear<T>> : std::true_type {};

template<typename T>
struct IsSTG : std::false_type {};

template<typename T>
struct IsSTG<property_types::ao_integrals::STG<T>> : std::true_type {};

template<typename T>
struct IsYukawa : std::false_type {};

template<typename T>
struct IsYukawa<property_types::ao_integrals::Yukawa<T>> : std::true_type {};

} // namespace detail_

template<typename T>
static constexpr auto op_v = detail_::LibintOp<T>::value;

template<typename T>
static constexpr auto is_doi_v = detail_::IsDOI<T>::value;

template<typename T>
static constexpr auto is_nuclear_v = detail_::IsNuclear<T>::value;

template<typename T>
static constexpr auto is_stg_v = detail_::IsSTG<T>::value;

template<typename T>
static constexpr auto is_yukawa_v = detail_::IsYukawa<T>::value;
} // namespace integrals
#pragma once
#include "integrals/property_types.hpp"

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

} // namespace detail_

template<typename T>
static constexpr auto op_v = detail_::LibintOp<T>::value;

} // namespace integrals
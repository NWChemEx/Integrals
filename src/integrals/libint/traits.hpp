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
struct NumberOfCenters;

template<typename T>
struct NumberOfCenters<pt::doi<T>> {
    static constexpr unsigned value = 4;
};

template<typename T>
struct IsDOI : std::false_type {};

template<typename T>
struct IsDOI<pt::doi<T>> : std::true_type {};

} // namespace detail_

template<typename T>
static constexpr auto op_v = detail_::LibintOp<T>::value;

template<typename T>
static constexpr auto number_of_centers_v = detail_::NumberOfCenters<T>::value;

template<typename T>
static constexpr auto is_doi_v = detail_::IsDOI<T>::value;

} // namespace integrals
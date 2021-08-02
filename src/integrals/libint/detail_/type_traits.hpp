#pragma once
#include <simde/types.hpp>

namespace integrals {
namespace detail_ {

template<typename T>
struct LibintOp;

// template<typename T>
// struct LibintOp<pt::doi<T>> {
//     static constexpr auto value = libint2::Operator::delta;
// };

template<>
struct LibintOp<simde::type::el_el_coulomb> {
    static constexpr auto value = libint2::Operator::coulomb;
};

template<>
struct LibintOp<simde::type::el_kinetic> {
    static constexpr auto value = libint2::Operator::kinetic;
};

template<>
struct LibintOp<simde::type::el_identity> {
    static constexpr auto value = libint2::Operator::overlap;
};

template<>
struct LibintOp<simde::type::el_nuc_coulomb> {
    static constexpr auto value = libint2::Operator::nuclear;
};

template<>
struct LibintOp<simde::type::el_el_stg> {
    static constexpr auto value = libint2::Operator::stg;
};

template<>
struct LibintOp<simde::type::el_el_yukawa> {
    static constexpr auto value = libint2::Operator::yukawa;
};

} // namespace detail_

template<typename T>
static constexpr auto op_v = detail_::LibintOp<T>::value;

} // namespace integrals

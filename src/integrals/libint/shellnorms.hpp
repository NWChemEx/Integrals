#pragma once
#include <libint2.hpp>
#include <pluginplay/module_base.hpp>

namespace integrals {

template<typename OperatorType>
DECLARE_MODULE(ShellNorms);

using ShellNormCoulomb = ShellNorms<simde::type::el_el_coulomb>;
using ShellNormSTG     = ShellNorms<simde::type::el_el_stg>;
using ShellNormYukawa  = ShellNorms<simde::type::el_el_yukawa>;

extern template class ShellNorms<simde::type::el_el_coulomb>;
extern template class ShellNorms<simde::type::el_el_stg>;
extern template class ShellNorms<simde::type::el_el_yukawa>;

} // namespace integrals

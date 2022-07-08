#pragma once
#include <pluginplay/module_base.hpp>
#include <simde/types.hpp>

namespace integrals {

template<std::size_t N, typename OperatorType>
DECLARE_MODULE(Libint);

extern template class Libint<2, simde::type::el_el_coulomb>;
extern template class Libint<3, simde::type::el_el_coulomb>;
extern template class Libint<4, simde::type::el_el_coulomb>;
extern template class Libint<2, simde::type::el_kinetic>;
extern template class Libint<2, simde::type::el_nuc_coulomb>;
extern template class Libint<2, simde::type::el_identity>;
extern template class Libint<2, simde::type::el_el_stg>;
extern template class Libint<3, simde::type::el_el_stg>;
extern template class Libint<4, simde::type::el_el_stg>;
extern template class Libint<2, simde::type::el_el_yukawa>;
extern template class Libint<3, simde::type::el_el_yukawa>;
extern template class Libint<4, simde::type::el_el_yukawa>;
extern template class Libint<2, simde::type::el_el_f12_commutator>;
extern template class Libint<3, simde::type::el_el_f12_commutator>;
extern template class Libint<4, simde::type::el_el_f12_commutator>;

DECLARE_MODULE(LibintDOI);

template<std::size_t L, typename OperatorType>
DECLARE_MODULE(LibintMultipole);

extern template class LibintMultipole<0, simde::type::el_dipole>;
extern template class LibintMultipole<1, simde::type::el_quadrupole>;
extern template class LibintMultipole<2, simde::type::el_octupole>;

} // namespace integrals

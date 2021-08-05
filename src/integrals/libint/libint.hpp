#pragma once
#include <pluginplay/module_base.hpp>
#include <simde/types.hpp>

namespace integrals {

template<std::size_t N, typename OperatorType>
DECLARE_MODULE(Libint);

DECLARE_MODULE(LibintDipole);
DECLARE_MODULE(LibintQuadrupole);
DECLARE_MODULE(LibintOctupole);
DECLARE_MODULE(LibintDOI);

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

} // namespace integrals

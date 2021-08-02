#pragma once
#include <pluginplay/module_base.hpp>
#include <simde/types.hpp>

namespace integrals {

template<std::size_t N, typename OperatorType>
DECLARE_MODULE(Libint);

// /// Multipoles require a fairly different backend than the other integral
// types
// /// so we specialize the Libint module for it.
// ///@{
// template<typename T>
// DECLARE_MODULE(Libint<pt::edipole<T>>);

// template<typename T>
// DECLARE_MODULE(Libint<pt::equadrupole<T>>);

// template<typename T>
// DECLARE_MODULE(Libint<pt::eoctopole<T>>);
// ///@}

// extern template class Libint<pt::doi<double>>;
// extern template class Libint<pt::edipole<double>>;
// extern template class Libint<pt::equadrupole<double>>;
// extern template class Libint<pt::eoctopole<double>>;
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

} // namespace integrals

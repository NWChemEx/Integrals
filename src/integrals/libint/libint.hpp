#pragma once
#include "integrals/property_types.hpp"
#include <sde/module_base.hpp>

namespace integrals {

template<typename PropType>
DECLARE_MODULE(Libint);

/// Multipoles require a fairly different backend than the other integral types
/// so we specialize the Libint module for it.
///@{
template<typename T>
DECLARE_MODULE(Libint<pt::edipole<T>>);

template<typename T>
DECLARE_MODULE(Libint<pt::equadrupole<T>>);

template<typename T>
DECLARE_MODULE(Libint<pt::eoctopole<T>>);
///@}

/// Typedefs of the modules
///@{
template<typename T>
using LibintDOI = Libint<pt::doi<T>>;

template<typename T>
using LibintEDipole = Libint<pt::edipole<T>>;

template<typename T>
using LibintEQuadrupole = Libint<pt::equadrupole<T>>;

template<typename T>
using LibintEOctopole = Libint<pt::eoctopole<T>>;

template<typename T>
using LibintERI2C = Libint<pt::eri2c<T>>;

template<typename T>
using LibintERI3C = Libint<pt::eri3c<T>>;

template<typename T>
using LibintERI4C = Libint<pt::eri4c<T>>;

template<typename T>
using LibintKinetic = Libint<pt::kinetic<T>>;

template<typename T>
using LibintNuclear = Libint<pt::nuclear<T>>;

template<typename T>
using LibintOverlap = Libint<pt::overlap<T>>;

template<typename T>
using LibintSTG2C = Libint<pt::stg2c<T>>;

template<typename T>
using LibintSTG3C = Libint<pt::stg3c<T>>;

template<typename T>
using LibintSTG4C = Libint<pt::stg4c<T>>;

template<typename T>
using LibintYukawa2C = Libint<pt::yukawa2c<T>>;

template<typename T>
using LibintYukawa3C = Libint<pt::yukawa3c<T>>;

template<typename T>
using LibintYukawa4C = Libint<pt::yukawa4c<T>>;
///@}

extern template class Libint<pt::doi<double>>;
extern template class Libint<pt::edipole<double>>;
extern template class Libint<pt::equadrupole<double>>;
extern template class Libint<pt::eoctopole<double>>;
extern template class Libint<pt::eri2c<double>>;
extern template class Libint<pt::eri3c<double>>;
extern template class Libint<pt::eri4c<double>>;
extern template class Libint<pt::kinetic<double>>;
extern template class Libint<pt::nuclear<double>>;
extern template class Libint<pt::overlap<double>>;
extern template class Libint<pt::stg2c<double>>;
extern template class Libint<pt::stg3c<double>>;
extern template class Libint<pt::stg4c<double>>;
extern template class Libint<pt::yukawa2c<double>>;
extern template class Libint<pt::yukawa3c<double>>;
extern template class Libint<pt::yukawa4c<double>>;

} // namespace integrals
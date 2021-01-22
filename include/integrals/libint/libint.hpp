#pragma once
#include "integrals/property_types.hpp"
#include <sde/module_base.hpp>

namespace integrals {

template<typename PropType>
DECLARE_MODULE(Libint);

template<typename T>
DECLARE_MODULE(Libint<pt::edipole<T>>);

template<typename T>
DECLARE_MODULE(Libint<pt::equadrupole<T>>);

template<typename T>
DECLARE_MODULE(Libint<pt::eoctopole<T>>);

template<typename T>
using LibintDOI = Libint<pt::doi<T>>;

template<typename T>
using LibintEDipole = Libint<pt::edipole<T>>;

template<typename T>
using LibintEQuadrupole = Libint<pt::equadrupole<T>>;

template<typename T>
using LibintEOctopole = Libint<pt::eoctopole<T>>;

extern template class Libint<pt::doi<double>>;
extern template class Libint<pt::edipole<double>>;
extern template class Libint<pt::equadrupole<double>>;
extern template class Libint<pt::eoctopole<double>>;

} // namespace integrals
#pragma once
#include <sde/module_base.hpp>

namespace integrals {

template<typename element_type = double>
DECLARE_MODULE(EDipoleInt);

template<typename element_type = double>
DECLARE_MODULE(EQuadrupoleInt);

template<typename element_type = double>
DECLARE_MODULE(EOctopoleInt);

extern template class EDipoleInt<double>;
extern template class EQuadrupoleInt<double>;
extern template class EOctopoleInt<double>;

using EDipole     = EDipoleInt<double>;
using EQuadrupole = EQuadrupoleInt<double>;
using EOctopole   = EOctopoleInt<double>;
} // namespace integrals
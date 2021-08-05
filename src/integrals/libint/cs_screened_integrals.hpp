#pragma once
#include <pluginplay/module_base.hpp>
#include <simde/tensor_representation/tensor_representation.hpp>

namespace integrals {

template<typename PropType>
DECLARE_MODULE(CauchySchwarzScreened);

extern template class CauchySchwarzScreened<simde::ERI3>;
extern template class CauchySchwarzScreened<simde::ERI4>;
extern template class CauchySchwarzScreened<simde::STG3>;
extern template class CauchySchwarzScreened<simde::STG4>;
extern template class CauchySchwarzScreened<simde::Yukawa3>;
extern template class CauchySchwarzScreened<simde::Yukawa4>;

} // namespace integrals

#pragma once
#include <pluginplay/module_base.hpp>
#include <simde/types.hpp>

namespace integrals {

template<std::size_t N, typename OperatorType>
DECLARE_MODULE(CSLibint);

extern template class CSLibint<3, simde::type::el_el_coulomb>;
extern template class CSLibint<4, simde::type::el_el_coulomb>;
extern template class CSLibint<3, simde::type::el_el_stg>;
extern template class CSLibint<4, simde::type::el_el_stg>;
extern template class CSLibint<3, simde::type::el_el_yukawa>;
extern template class CSLibint<4, simde::type::el_el_yukawa>;

} // namespace integrals

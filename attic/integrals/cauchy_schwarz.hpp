#pragma once
#include <libint2.hpp>
#include <sde/module_base.hpp>

namespace integrals {

template<typename element_type = double,
         libint2::Operator op  = libint2::Operator::coulomb>
DECLARE_MODULE(CauchySchwarz);

extern template class CauchySchwarz<double>;
extern template class CauchySchwarz<double, libint2::Operator::stg>;
extern template class CauchySchwarz<double, libint2::Operator::yukawa>;

using CS_ERI    = CauchySchwarz<double>;
using CS_STG    = CauchySchwarz<double, libint2::Operator::stg>;
using CS_Yukawa = CauchySchwarz<double, libint2::Operator::yukawa>;

} // namespace integrals
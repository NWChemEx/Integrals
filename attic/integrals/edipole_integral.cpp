#include "integrals/emultipole_integrals.hpp"
#include "integrals/libint_integral.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include <chemist/ta_helpers/einsum/einsum.hpp>
#include <property_types/ao_integrals/emultipole.hpp>
#include <property_types/ao_integrals/overlap.hpp>

namespace integrals {

template<typename element_type>
using overlap_type = property_types::ao_integrals::Overlap<element_type>;
template<typename element_type>
using eDipole_type = property_types::ao_integrals::EDipole<element_type>;
template<typename element_type>
using libint_type = property_types::LibIntIntegral<element_type>;
template<typename element_type>
using tensor = typename type::tensor<element_type>;
template<typename element_type>
using value_type = typename tensor<element_type>::value_type;

template<typename element_type>
EDipoleInt<element_type>::EDipoleInt() : sde::ModuleBase(this) {
    description("Computes dipole integrals with Libint");
    satisfies_property_type<overlap_type<element_type>>();
    satisfies_property_type<eDipole_type<element_type>>();
    satisfies_property_type<libint_type<element_type>>();
}

template<typename element_type>
sde::type::result_map EDipoleInt<element_type>::run_(
  sde::type::input_map inputs, sde::type::submodule_map submods) const {}

template class EDipoleInt<double>;

} // namespace integrals

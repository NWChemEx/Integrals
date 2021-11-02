#include <simde/simde.hpp>

namespace integrals {

template<std::size_t N, typename OpType>
DECLARE_MODULE(StandardTransform);

extern template class StandardTransform<2, simde::type::el_scf_k>;
extern template class StandardTransform<2, simde::type::fock>;
extern template class StandardTransform<3, simde::type::el_el_coulomb>;
extern template class StandardTransform<4, simde::type::el_el_coulomb>;
extern template class StandardTransform<4, simde::type::el_el_f12_commutator>;
extern template class StandardTransform<4, simde::type::el_el_stg>;
extern template class StandardTransform<4, simde::type::el_el_yukawa>;
} // namespace integrals

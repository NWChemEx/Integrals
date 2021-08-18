#include <simde/simde.hpp>

namespace integrals {

template<std::size_t N, typename OpType>
DECLARE_MODULE(StandardTransform);

extern template class StandardTransform<4, simde::type::el_el_coulomb>;

} // namespace integrals

#include "integrals/transformed.hpp"

namespace integrals {

template class Transformed<pt::doi<double>>;
template class Transformed<pt::correlation_factor_4c<double>>;
template class Transformed<pt::eri2c<double>>;

} // namespace integrals
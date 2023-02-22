#include "libint_factory.hpp"

namespace integrals::libint {

template<typename std::size_t N, typename OperatorType>
TEMPLATED_MODULE_RUN(MakeLibintFactory, N, OperatorType) {
    /// Grab input information
    auto bases  = unpack_bases<N>(inputs);
    auto op_str = OperatorType().as_string();
    auto op     = inputs.at(op_str).template value<const OperatorType&>();
    auto thresh = inputs.at("Threshold").value<double>();
    auto deriv  = 0; // TODO: Get from TMP

    auto pfactory =
      std::make_unique<LibintFactory>(std::move(bases), op, thresh, deriv));
    IntegralFactory fac(std::move(pfactory));
}

} // namespace integrals::libint
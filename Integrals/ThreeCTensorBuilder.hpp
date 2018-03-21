#pragma once
#include <Integrals/IntegralTensorBuilder.hpp>
#include <Integrals/nwx_libint/LibInt3C.hpp>

namespace Integrals {

template<typename integral_type>
class ThreeCTensorBuilder : public IntegralTensorBuilder<3> {
public:
    using TensorType = IntegralTensorBuilder<3>::TensorType;

    ThreeCTensorBuilder() = default;
    ~ThreeCTensorBuilder() = default;

    std::vector<TensorType> compute(const LibChemist::Molecule &molecule,
             const std::vector<LibChemist::BasisSet> &basissets) const override;
};

extern template class ThreeCTensorBuilder<nwx_libint::DF3C2E>;

}

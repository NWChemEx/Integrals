#pragma once
#include <IntegralsEx/IntegralTensorBuilder.hpp>
#include <IntegralsEx/nwx_libint/LibInt3C.hpp>

namespace IntegralsEx {

template<typename integral_type>
class ThreeCTensorBuilder : public IntegralTensorBuilder {
public:
    using TensorType = IntegralTensorBuilder::TensorType;

    ThreeCTensorBuilder() = default;
    ~ThreeCTensorBuilder() = default;

    std::vector<TensorType> compute(const LibChemist::SetOfAtoms &atoms,
             const std::vector<LibChemist::BasisSet> &basissets) const override;
};

extern template class ThreeCTensorBuilder<nwx_libint::DF3C2E>;

}

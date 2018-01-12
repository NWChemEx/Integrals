#pragma once
#include <IntegralsEx/IntegralTensorBuilder.hpp>
#include <IntegralsEx/nwx_libint/LibInt4C.hpp>

namespace IntegralsEx {

template<typename integral_type>
class FourCTensorBuilder : public IntegralTensorBuilder<4> {
public:
    using TensorType = IntegralTensorBuilder<4>::TensorType;

    FourCTensorBuilder() = default;
    ~FourCTensorBuilder() = default;

    std::vector<TensorType> compute(const LibChemist::SetOfAtoms &atoms,
             const std::vector<LibChemist::BasisSet> &basissets) const override;
};

extern template class FourCTensorBuilder<nwx_libint::ERI>;
extern template class FourCTensorBuilder<nwx_libint::Delta>;

}

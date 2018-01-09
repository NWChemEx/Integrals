#pragma once
#include <IntegralsEx/IntegralTensorBuilder.hpp>
#include <IntegralsEx/nwx_libint/LibInt2C.hpp>

namespace IntegralsEx {

template<typename integral_type>
class TwoCTensorBuilder : public IntegralTensorBuilder<2> {
public:
    using TensorType = IntegralTensorBuilder<2>::TensorType;

    TwoCTensorBuilder() = default;
    ~TwoCTensorBuilder() = default;

    std::vector<TensorType> compute(const LibChemist::SetOfAtoms &atoms,
             const std::vector<LibChemist::BasisSet> &basissets) const override;
};

extern template class TwoCTensorBuilder<nwx_libint::Overlap>;
extern template class TwoCTensorBuilder<nwx_libint::Kinetic>;
extern template class TwoCTensorBuilder<nwx_libint::NuclearElectron>;
extern template class TwoCTensorBuilder<nwx_libint::Metric>;

}

#pragma once
#include <Integrals/IntegralTensorBuilder.hpp>
#include <Integrals/nwx_libint/LibInt4C.hpp>

namespace Integrals {

template<typename integral_type>
class FourCTensorBuilder : public IntegralTensorBuilder<4> {
public:
    using TensorType = IntegralTensorBuilder<4>::TensorType;

    FourCTensorBuilder() = default;
    ~FourCTensorBuilder();
  

    std::vector<TensorType> compute(const LibChemist::Molecule &molecule,
             const std::vector<LibChemist::BasisSet> &basissets) const override;
};

  //External destructor definition forces compiler to produce the destructor for extern template classes
  template<typename integral_type>
  FourCTensorBuilder<integral_type>::~FourCTensorBuilder() = default;

extern template class FourCTensorBuilder<nwx_libint::ERI>;
extern template class FourCTensorBuilder<nwx_libint::Delta>;

}

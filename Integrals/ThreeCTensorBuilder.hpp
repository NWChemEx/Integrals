#pragma once
#include <Integrals/IntegralTensorBuilder.hpp>
#include <Integrals/nwx_libint/LibInt3C.hpp>

namespace Integrals {

template<typename integral_type>
class ThreeCTensorBuilder : public IntegralTensorBuilder<3> {
public:
    using TensorType = IntegralTensorBuilder<3>::TensorType;

    ThreeCTensorBuilder() = default;
    ~ThreeCTensorBuilder();

    std::vector<TensorType> compute(const LibChemist::Molecule &molecule,
             const std::vector<LibChemist::BasisSet> &basissets) const override;
};

  //External destructor definition forces compiler to produce the destructor for extern template classes
  template<typename integral_type>
  ThreeCTensorBuilder<integral_type>::~ThreeCTensorBuilder() = default;

extern template class ThreeCTensorBuilder<nwx_libint::DF3C2E>;

}

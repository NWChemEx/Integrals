#pragma once
#include <Integrals/IntegralTensorBuilder.hpp>
#include <Integrals/nwx_libint/LibInt2C.hpp>

namespace Integrals {

template<typename integral_type>
class TwoCTensorBuilder : public IntegralTensorBuilder<2> {
public:
    using TensorType = IntegralTensorBuilder<2>::TensorType;

    TwoCTensorBuilder() = default;
    ~TwoCTensorBuilder();

    std::vector<TensorType> compute(const LibChemist::Molecule &molecule,
             const std::vector<LibChemist::BasisSet> &basissets) const override;
};

  //External destructor definition forces compiler to produce the destructor for extern template classes
  template<typename integral_type>
  TwoCTensorBuilder<integral_type>::~TwoCTensorBuilder() = default;


extern template class TwoCTensorBuilder<nwx_libint::Overlap>;
extern template class TwoCTensorBuilder<nwx_libint::Kinetic>;
extern template class TwoCTensorBuilder<nwx_libint::NuclearElectron>;
extern template class TwoCTensorBuilder<nwx_libint::Metric>;
extern template class TwoCTensorBuilder<nwx_libint::EDipole>;
extern template class TwoCTensorBuilder<nwx_libint::EQuadrupole>;
extern template class TwoCTensorBuilder<nwx_libint::EOctopole>;

}

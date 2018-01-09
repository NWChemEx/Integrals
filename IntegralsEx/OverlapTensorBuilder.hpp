#pragma once
#include <IntegralsEx/IntegralTensorBuilder.hpp>

namespace IntegralsEx {

class OverlapTensorBuilder : public IntegralTensorBuilder {
public:
    OverlapTensorBuilder() = default;
    ~OverlapTensorBuilder() = default;

    std::vector<TensorType> compute(const LibChemist::SetOfAtoms &atoms,
             const std::vector<LibChemist::BasisSet> &basissets) const override;
};

}

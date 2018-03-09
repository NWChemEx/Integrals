#pragma once
#include <LibChemist/SetOfAtoms.hpp>
#include <unsupported/Eigen/CXX11/Tensor>

namespace Integrals{

template<size_t order>
class IntegralTensorBuilder
{
public:
     using TensorType = Eigen::Tensor<double,order>;

     IntegralTensorBuilder() = default;
     virtual ~IntegralTensorBuilder() {}

     virtual std::vector<TensorType> compute(const LibChemist::SetOfAtoms &atoms,
             const std::vector<LibChemist::BasisSet> &basissets) const = 0;
};

}

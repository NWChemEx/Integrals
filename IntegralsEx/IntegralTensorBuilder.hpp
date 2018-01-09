#pragma once
#include <LibChemist/SetOfAtoms.hpp>
#include <unsupported/Eigen/CXX11/Tensor>
namespace IntegralsEx{


class IntegralTensorBuilder
{
public:
     using TensorType = Eigen::Tensor<double,2>;

     IntegralTensorBuilder() = default;
     virtual ~IntegralTensorBuilder() = default;

     virtual std::vector<TensorType> compute(const LibChemist::SetOfAtoms &atoms,
             const std::vector<LibChemist::BasisSet> &basissets) const = 0;
};

}

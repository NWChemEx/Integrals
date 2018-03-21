#pragma once
#include <LibChemist/Molecule.hpp>
#include <unsupported/Eigen/CXX11/Tensor>

namespace Integrals{

template<size_t order>
class IntegralTensorBuilder
{
public:
     using TensorType = Eigen::Tensor<double,order>;

     IntegralTensorBuilder() = default;
     virtual ~IntegralTensorBuilder() {}

     virtual std::vector<TensorType> compute(const LibChemist::Molecule &molecule,
             const std::vector<LibChemist::BasisSet> &basissets) const = 0;
};

}

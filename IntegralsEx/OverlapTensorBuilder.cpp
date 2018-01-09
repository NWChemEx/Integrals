#include <stdexcept>
#include <IntegralsEx/OverlapTensorBuilder.hpp>
#include <IntegralsEx/nwx_libint/LibInt2C.hpp>

namespace IntegralsEx {

std::vector<typename OverlapTensorBuilder::TensorType> compute(const LibChemist::SetOfAtoms &atoms,
             const std::vector<LibChemist::BasisSet> &basissets) const {
    if (basissets.size() != 2)
        throw std::length_error("Wrong number of basis sets");

    nwx_libint::Overlap libints(0,atoms,basissets[0],basissets[1]);
    IntegralsEx::TwoCenterIntegral *Ints = &libints;
    OverlapTensorBuilder::TensorType tensor({basissets[0].size(),basissets[1].size()});
    tensor.setZero();

    std::size_t off_i = 0;
    for(std::size_t i = 0; i < basissets[0].ngens.size(); i++) {
        const std::size_t ni = basisets[0].shellsize(i);
        std::size_t off_j = 0;
        for(std::size_t j = 0; j < basissets[1].ngens.size(); j++) {
            const std::size_t nj = basisets[1].shellsize(j);
            const double* buffer=Ints->calculate(i,j);
            if(buffer==nullptr)continue;
            for (std::size_t si = 0; si < ni; si++) {
                for (std::size_t sj = 0; sj < nj; sj++) {
                    tensor(off_i + si, off_j + sj) = buffer++;
                }
            }
            off_j += nj;
        }
        off_i += ni;
    }
    return {std::move(tensor)};
}

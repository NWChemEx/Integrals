#include <stdexcept>
#include <IntegralsEx/ThreeCTensorBuilder.hpp>

namespace IntegralsEx {

template<typename libint_type> 
std::vector<typename IntegralTensorBuilder<3>::TensorType> kernel(const LibChemist::SetOfAtoms &atoms,
             const std::vector<LibChemist::BasisSet> &basissets) {
    if (basissets.size() != 3)
        throw std::length_error("Wrong number of basis sets");

    libint_type libints(0,atoms,basissets[0],basissets[1],basissets[2]);
    IntegralsEx::ThreeCenterIntegral *Ints = &libints;
    IntegralTensorBuilder<3>::TensorType rv(std::array<long int,3>{basissets[0].size(),basissets[1].size(),basissets[2].size()});
    rv.setZero();

    std::size_t off_i = 0;
    for(std::size_t i = 0; i < basissets[0].ngens.size(); i++) {
        const std::size_t ni = basissets[0].shellsize(i);
        std::size_t off_j = 0;
        for(std::size_t j = 0; j < basissets[1].ngens.size(); j++) {
            const std::size_t nj = basissets[1].shellsize(j);
            std::size_t off_k = 0;
            for(std::size_t k = 0; k < basissets[2].ngens.size(); k++) {
                const std::size_t nk = basissets[2].shellsize(k);
                const double* buffer=Ints->calculate(i,j,k);
                if(buffer==nullptr)continue;
                for (std::size_t si = 0; si < ni; si++) {
                    for (std::size_t sj = 0; sj < nj; sj++) {
                        for (std::size_t sk = 0; sk < nk; sk++) {
                            rv(off_i + si, off_j + sj, off_k + sk) = *buffer++;
                        }
                    }
                }
                off_k += nk;
            }
            off_j += nj;
        }
        off_i += ni;
    }
    return {std::move(rv)};
}

template<typename libint_type>
std::vector<typename IntegralTensorBuilder<3>::TensorType> ThreeCTensorBuilder<libint_type>::compute(const LibChemist::SetOfAtoms &atoms,
             const std::vector<LibChemist::BasisSet> &basissets) const {
    return kernel<libint_type>(atoms, basissets);
}

template class ThreeCTensorBuilder<nwx_libint::DF3C2E>;

}

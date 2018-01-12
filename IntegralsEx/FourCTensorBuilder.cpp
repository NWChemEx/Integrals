#include <stdexcept>
#include <IntegralsEx/FourCTensorBuilder.hpp>

namespace IntegralsEx {

template<typename libint_type> 
std::vector<typename IntegralTensorBuilder<4>::TensorType> kernel(const LibChemist::SetOfAtoms &atoms,
             const std::vector<LibChemist::BasisSet> &basissets) {
    if (basissets.size() != 4)
        throw std::length_error("Wrong number of basis sets");

    libint_type libints(0,atoms,basissets[0],basissets[1],basissets[2],basissets[3]);
    IntegralsEx::FourCenterIntegral *Ints = &libints;
    IntegralTensorBuilder<4>::TensorType rv(std::array<long int,4>{basissets[0].size(),basissets[1].size(),basissets[2].size(),basissets[3].size()});
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
                std::size_t off_l = 0;
                for(std::size_t l = 0; l < basissets[3].ngens.size(); l++) {
                    const std::size_t nl = basissets[3].shellsize(l);
                    const double* buffer=Ints->calculate(i,j,k,l);
                    if(buffer==nullptr)continue;
                    for (std::size_t si = 0; si < ni; si++) {
                        for (std::size_t sj = 0; sj < nj; sj++) {
                            for (std::size_t sk = 0; sk < nk; sk++) {
                                for (std::size_t sl = 0; sl < nl; sl++) 
                                    rv(off_i + si, off_j + sj, off_k + sk, off_l +sl) = *buffer++;
                            }   
                        }
                    }
                    off_l += nl;
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
std::vector<typename IntegralTensorBuilder<4>::TensorType> FourCTensorBuilder<libint_type>::compute(const LibChemist::SetOfAtoms &atoms,
             const std::vector<LibChemist::BasisSet> &basissets) const {
    return kernel<libint_type>(atoms, basissets);
}

template class FourCTensorBuilder<nwx_libint::ERI>;
template class FourCTensorBuilder<nwx_libint::Delta>;

}

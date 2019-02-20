#include <stdexcept>
#include <Integrals/FourCTensorBuilder.hpp>

namespace Integrals {

template<typename libint_type> 
std::vector<typename IntegralTensorBuilder<4>::TensorType> kernel(const LibChemist::Molecule &molecule,
             const std::vector<LibChemist::BasisSet> &basissets) {
    if (basissets.size() != 4)
        throw std::length_error("Wrong number of basis sets");

    libint_type libints(0,molecule,basissets[0],basissets[1],basissets[2],basissets[3]);
    Integrals::FourCenterIntegral *Ints = &libints;
    IntegralTensorBuilder<4>::TensorType init_tensor(std::array<long int,4>
                                                     {static_cast<long>(basissets[0].size()),static_cast<long>(basissets[1].size()),
                                                      static_cast<long>(basissets[2].size()),static_cast<long>(basissets[3].size())});
    init_tensor.setZero();
    std::vector<IntegralTensorBuilder<4>::TensorType> rv(Ints->n_components(),init_tensor);

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
                    std::vector<const double*> buf_vec=Ints->calculate(i,j,k,l);
                    size_t counter = 0;
                    for (auto buffer : buf_vec) {
                        if(buffer==nullptr)continue;
                        for (std::size_t si = 0; si < ni; si++) {
                            for (std::size_t sj = 0; sj < nj; sj++) {
                                for (std::size_t sk = 0; sk < nk; sk++) {
                                    for (std::size_t sl = 0; sl < nl; sl++) 
                                        rv[counter](off_i + si, off_j + sj, 
                                                    off_k + sk, off_l +sl) = *buffer++;
                                }   
                            }
                        }
                    counter++;
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
std::vector<typename IntegralTensorBuilder<4>::TensorType> FourCTensorBuilder<libint_type>::compute(const LibChemist::Molecule &molecule,
             const std::vector<LibChemist::BasisSet> &basissets) const {
    return kernel<libint_type>(molecule, basissets);
}

template class FourCTensorBuilder<nwx_libint::ERI>;
template class FourCTensorBuilder<nwx_libint::Delta>;

}

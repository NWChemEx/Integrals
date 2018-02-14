#include <stdexcept>
#include <IntegralsEx/TwoCTensorBuilder.hpp>

namespace IntegralsEx {

template<typename libint_type> 
std::vector<typename IntegralTensorBuilder<2>::TensorType> kernel(const LibChemist::SetOfAtoms &atoms,
             const std::vector<LibChemist::BasisSet> &basissets) {
    if (basissets.size() != 2)
        throw std::length_error("Wrong number of basis sets");

    libint_type libints(0,atoms,basissets[0],basissets[1]);
    IntegralsEx::TwoCenterIntegral *Ints = &libints;
    IntegralTensorBuilder<2>::TensorType init_tensor(std::array<long int,2>
                                                    {basissets[0].size(),basissets[1].size()});
    init_tensor.setZero();
    std::vector<IntegralTensorBuilder<2>::TensorType> rv(Ints->n_components(),init_tensor);

    std::size_t off_i = 0;
    for(std::size_t i = 0; i < basissets[0].ngens.size(); i++) {
        const std::size_t ni = basissets[0].shellsize(i);
        std::size_t off_j = 0;
        for(std::size_t j = 0; j < basissets[1].ngens.size(); j++) {
            const std::size_t nj = basissets[1].shellsize(j);
            std::vector<const double*> buf_vec=Ints->calculate(i,j);
            size_t counter = 0;
            for (auto buffer : buf_vec) {
                if(buffer==nullptr)continue;
                for (std::size_t si = 0; si < ni; si++) {
                    for (std::size_t sj = 0; sj < nj; sj++) {
                        rv[counter](off_i + si, off_j + sj) = *buffer++;
                    }
                }
                counter++;
            }
            off_j += nj;
        }
        off_i += ni;
    }
    return rv;
}

template<typename libint_type>
std::vector<typename IntegralTensorBuilder<2>::TensorType> TwoCTensorBuilder<libint_type>::compute(const LibChemist::SetOfAtoms &atoms,
             const std::vector<LibChemist::BasisSet> &basissets) const {
    return kernel<libint_type>(atoms, basissets);
}

template class TwoCTensorBuilder<nwx_libint::Overlap>;
template class TwoCTensorBuilder<nwx_libint::Kinetic>;
template class TwoCTensorBuilder<nwx_libint::NuclearElectron>;
template class TwoCTensorBuilder<nwx_libint::Metric>;
template class TwoCTensorBuilder<nwx_libint::EDipole>;
template class TwoCTensorBuilder<nwx_libint::EQuadrupole>;
template class TwoCTensorBuilder<nwx_libint::EOctopole>;

}

#include "Integrals/nwx_libint/nwx_libint.hpp"

namespace nwx_libint{

libint2::BasisSet make_basis(const LibChemist::AOBasisSet& bs)
{
    using Contraction=libint2::Shell::Contraction;
    libint2::BasisSet shells;
    for(const auto& shelli : bs) {
        const auto nprims = shelli.nprims();
        std::vector<double> alphas(&shelli.alpha(0), &shelli.alpha(nprims));
        std::vector<double> coefs(&shelli.coef(0), &shelli.coef(nprims));
        Contraction cont({shelli.l(), shelli.pure(), coefs});
        std::vector<Contraction> conts{cont};
        shells.push_back(libint2::Shell(alphas, conts, shelli.center()));
    }
    return shells;
}

}//End namespace

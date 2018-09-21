#include "Integrals/nwx_libint/nwx_libint.hpp"

namespace nwx_libint{

libint2::BasisSet make_basis(const LibChemist::AOBasisSet& bs)
{
    using Contraction=libint2::Shell::Contraction;
    libint2::BasisSet shells;
    for(const auto& shelli : bs) {
        const auto nprims = shelli.nprims();
        libint2::svector<double> alphas(&shelli.alpha(0), &shelli.alpha(nprims));
        libint2::svector<double> coefs(&shelli.coef(0), &shelli.coef(nprims));
        int l = shelli.l();
        Contraction cont({l, shelli.pure(), coefs});
        libint2::svector<Contraction> conts{cont};
        shells.push_back(libint2::Shell(alphas, conts, shelli.center()));
    }
    return shells;
}

}//End namespace

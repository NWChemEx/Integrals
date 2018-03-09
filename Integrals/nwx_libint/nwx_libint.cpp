#include "Integrals/nwx_libint/nwx_libint.hpp"

namespace nwx_libint{

libint2::BasisSet make_basis(const LibChemist::BasisSet& bs)
{
    using Contraction=libint2::Shell::Contraction;
    const auto pure=LibChemist::ShellType::SphericalGaussian;
    libint2::BasisSet shells;
    size_t j = 0;
    auto alpha_start = bs.alphas.begin();
    auto coef_start = bs.coefs.begin();
    for(size_t i=0; i<bs.nprims.size(); ++i)
    {
        const double x=bs.centers[j++],
                     y=bs.centers[j++],
                     z=bs.centers[j++];
        const std::array<double,3> carts({x,y,z});
        auto alpha_end = alpha_start + bs.nprims[i];
        //LibInt2 does not support general contractions
        for(size_t k=0; k<bs.ngens[i]; ++k)
        {
            int L;
            std::vector<Contraction> conts;
            if (bs.ls[i] < 0) {
                L=k;
            } else {
                L=bs.ls[i];
            }
            const bool is_pure=bs.types[i]==pure;
            auto coef_end = coef_start + bs.nprims[i];
            conts.push_back(Contraction({L,is_pure,std::vector<double>(coef_start,coef_end)}));
            coef_start = coef_end;
            shells.push_back(libint2::Shell(std::vector<double>(alpha_start,alpha_end),conts, carts));
        }
        alpha_start = alpha_end;
    }
    return shells;
}

}//End namespace

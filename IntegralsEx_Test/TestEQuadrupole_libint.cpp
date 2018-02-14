#include <IntegralsEx/nwx_libint/LibInt2C.hpp> 
#include "TestCommon.hpp"
#include "H2O_STO3G_Multipole.hpp"
#include <iostream>

TEST_CASE("Testing EQuadrupole class"){

    auto atoms=make_atoms();
    auto bs=get_basis("PRIMARY",atoms);
    nwx_libint::EQuadrupole libints(0,atoms,bs,bs);
    IntegralsEx::TwoCenterIntegral *Ints = &libints;

    size_t off_i = 0;
    for(size_t i=0;i<5;++i)
    {
        const size_t ni = bs.shellsize(i);
        size_t off_j = 0;
        for(size_t j=0;j<5;++j)
        {
            const size_t nj = bs.shellsize(j);
            std::vector<const double*> buf_vec=Ints->calculate(i,j);
            size_t ncomp = 0;
            for (auto buffer : buf_vec){
                if(buf_vec[0]==nullptr)continue;
                for (size_t si = 0; si < ni; si++) {
                    for (size_t sj = 0; sj < nj; sj++) {
                        REQUIRE(*buffer++ == Approx(corr[ncomp][7*(off_i+si)+(off_j+sj)]).epsilon(eps).margin(marg));
                    }
                }
                ncomp++;
            }
            off_j += nj;
        }
        off_i += ni;
    }
}

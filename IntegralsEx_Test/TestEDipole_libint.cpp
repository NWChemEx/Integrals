#include <IntegralsEx/nwx_libint/LibInt2C.hpp> 
#include "TestCommon.hpp"
#include "H2O_STO3G_Multipole.hpp"

TEST_CASE("Testing EDipole class"){

    auto atoms=make_atoms();
    auto bs=get_basis("PRIMARY",atoms);
    nwx_libint::EDipole libints(0,atoms,bs,bs);
    IntegralsEx::TwoCenterIntegral *Ints = &libints;

    size_t counter=0;
    for(size_t si=0;si<5;++si)
    {
        const size_t ni = bs.shellsize(si);
        for(size_t sj=0;sj<5;++sj)
        {
            const size_t nj = bs.shellsize(sj);
            std::vector<const double*> buf_vec=Ints->calculate(si,sj);
            for (size_t buf_idx = 0; buf_idx < ni*nj; buf_idx++) {
                if(buf_vec[0]==nullptr)continue;
            //    REQUIRE(buf_vec[0][buf_idx] == Approx(corr[0][counter]).epsilon(eps).margin(marg));
                counter++;
            }
        }
    }
}

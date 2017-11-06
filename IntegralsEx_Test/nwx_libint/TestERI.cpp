#include <IntegralsEx/nwx_libint/LibInt4C.hpp>
#include "../TestCommon.hpp"
#define CATCH_CONFIG_MAIN
#include "catch/catch.hpp"
#include "../H2O_STO3G_ERI.hpp"


//Computes the ERI integrals for water in STO-3G
TEST_CASE("Testing the 4C ERI class"){

    auto atoms=make_atoms();
    auto bs=get_basis("PRIMARY",atoms);
   
    nwx_libint::ERI libints(0,atoms,bs,bs,bs,bs);
    IntegralsEx::FourCenterIntegral *eri = &libints;

    size_t counter=0;
    for(size_t si=0;si<5;++si)
    {
        const size_t ni=bs.shellsize(si);
        for(size_t sj=0;sj<=si;++sj)
        {
            const size_t nj=bs.shellsize(sj);
            for(size_t sk=0;sk<=si;++sk)
            {
                const size_t nk=bs.shellsize(sk);
                for(size_t sl=0;sl<=(si==sk?sj:sk);++sl)
                {
                    const size_t nl=bs.shellsize(sl);
                    const double* buffer=eri->calculate(si,sj,sk,sl);
                    if(buffer==nullptr)continue;
                    ptr_wrapper wrapped_buffer={buffer,ni*nj*nk*nl};
                    REQUIRE(std::equal(corr[counter].begin(),corr[counter].end(),
                               wrapped_buffer.begin(),wrapped_buffer.end(),
                               [](double x, double y)->bool {return x==Approx(y);}));
                    ++counter;
                }
            }
        }
    }
}

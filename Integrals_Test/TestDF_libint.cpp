#include <Integrals/nwx_libint/LibInt3C.hpp>

#include "TestCommon.hpp"
#include "H2O_STO3G_DF.hpp"

//Computes the three-center, two-electron integrals for water in sto-3g
//Note: as with the metric tensor using sto-3g instead of a proper fitting basis
//should be fine for testing purposes.

TEST_CASE("Testing DF3C2E"){

    auto molecule=make_molecule();
    auto bs=molecule.get_basis("sto-3g");

    nwx_libint::DF3C2E libints(0,molecule,bs,bs,bs);
    Integrals::ThreeCenterIntegral *df3c = &libints;

    size_t counter=0;
    const size_t nshells=5;
    //TODO: sk<=sj
    for(size_t si=0;si<nshells;++si)
    {
        const size_t ni=bs.shellsize(si);
        for(size_t sj=0;sj<nshells;++sj)
        {
            const size_t nj=bs.shellsize(sj);
            for(size_t sk=0;sk<nshells;++sk)
            {
                    const size_t nk=bs.shellsize(sk);
                    const double* buffer=(df3c->calculate(si,sj,sk))[0];
                    if(buffer==nullptr)continue;
                    ptr_wrapper wrapped_buffer={buffer,ni*nj*nk};
                    compare_integrals(wrapped_buffer, corr[counter]);
                    ++counter;
            }
        }
    }
}

#include <IntegralsEx/nwx_libint/LibInt2C.hpp> 
#define CATCH_CONFIG_MAIN
#include "TestCommon.hpp"

//Computes the overlap integrals for water in STO-3G
std::vector<std::vector<double>> corr={
    {//(  0 |  0 )
    1.0000000000000004,},
    {//(  1 |  0 )
    0.2367039365108476,},
    {//(  1 |  1 )
    1.0000000000000002,},
    {//(  2 |  0 )
    0.0000000000000000,0.0000000000000000,0.0000000000000000,},
    {//(  2 |  1 )
    0.0000000000000000,-0.0000000000000000,0.0000000000000000,},
    {//(  2 |  2 )
    1.0000000000000002,0.0000000000000000,0.0000000000000000,
    0.0000000000000000,1.0000000000000002,0.0000000000000000,
    0.0000000000000000,0.0000000000000000,1.0000000000000002,},
    {//(  3 |  0 )
    0.0384055905135490,},
    {//(  3 |  1 )
    0.3861387813310925,},
    {//(  3 |  2 )
    0.2684376412681760,0.2097276494226497,0.0000000000000000,},
    {//(  3 |  3 )
    1.0000000000000002,},
    {//(  4 |  0 )
    0.0384055905135490,},
    {//(  4 |  1 )
    0.3861387813310925,},
    {//(  4 |  2 )
    -0.2684376412681760,0.2097276494226497,0.0000000000000000,},
    {//(  4 |  3 )
    0.1817608668218927,},
    {//(  4 |  4 )
    1.0000000000000002,},
};

TEST_CASE("Testing Overlap class"){

    auto atoms=make_atoms();
    auto bs=get_basis("PRIMARY",atoms);
    nwx_libint::Overlap libints(0,atoms,bs,bs);
    IntegralsEx::TwoCenterIntegral *Ints = &libints;

    size_t counter=0;
    for(size_t si=0;si<5;++si)
    {
        const size_t ni = bs.shellsize(si);
        for(size_t sj=0;sj<=si;++sj)
        {
            const size_t nj = bs.shellsize(sj);
            const double* buffer=Ints->calculate(si,sj);
            if(buffer==nullptr)continue;
            ptr_wrapper wrapped_buffer={buffer,ni*nj};
            compare_integrals(wrapped_buffer, corr[counter]);
            ++counter;
        }
    }
}


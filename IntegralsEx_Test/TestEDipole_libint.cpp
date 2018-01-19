#include <IntegralsEx/nwx_libint/LibInt2C.hpp> 

#include "TestCommon.hpp"
#include <iomanip>
#include <iostream>

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

TEST_CASE("Testing EDipole class"){

    auto atoms=make_atoms();
    auto bs=get_basis("PRIMARY",atoms);
    nwx_libint::EDipole libints(0,atoms,bs,bs);
    IntegralsEx::TwoCenterIntegral *Ints = &libints;

    size_t counter=0;
    for(size_t si=0;si<5;++si)
    {
        const size_t ni = bs.shellsize(si);
        for(size_t sj=0;sj<=si;++sj)
        {
            const size_t nj = bs.shellsize(sj);
            std::vector<const double*> buf_vec=Ints->calculate(si,sj);
            for (auto buffer : buf_vec) {
                if(buffer==nullptr)continue;
                size_t k = 0;
                for(size_t i=0;i<ni;++i)
                {
                    for(size_t j=0;j<nj;++j)
                        std::cout << std::fixed << std::setprecision(8) << buffer[k++] << ",";
                    std::cout << std::endl;
                }
                std::cout << std::endl;
                if (buffer == buf_vec[0]) {
                    ptr_wrapper wrapped_buffer={buffer,ni*nj};
                    compare_integrals(wrapped_buffer, corr[counter]);
                    ++counter;
                }
            }
        }
    }
}


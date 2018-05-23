#include <Integrals/LibintIntegral.hpp>

#include "TestCommon.hpp"

//Computes the kinetic energy integrals for water in STO-3G
std::vector<std::vector<double>> corr={
    {//(  0 |  0 )
    29.0031999455395848,},
    {//(  1 |  0 )
    -0.1680109393164923,},
    {//(  1 |  1 )
    0.8081279549303477,},
    {//(  2 |  0 )
    0.0000000000000000,0.0000000000000001,0.0000000000000000,},
    {//(  2 |  1 )
    0.0000000000000000,-0.0000000000000000,0.0000000000000000,},
    {//(  2 |  2 )
    2.5287311981947642,0.0000000000000000,0.0000000000000000,
    0.0000000000000000,2.5287311981947642,0.0000000000000000,
    0.0000000000000000,0.0000000000000000,2.5287311981947642,},
    {//(  3 |  0 )
    -0.0084163851854474,},
    {//(  3 |  1 )
    0.0705173385189986,},
    {//(  3 |  2 )
    0.1470905524127554,0.1149203802569079,0.0000000000000000,},
    {//(  3 |  3 )
    0.7600318835666090,},
    {//(  4 |  0 )
    -0.0084163851854474,},
    {//(  4 |  1 )
    0.0705173385189986,},
    {//(  4 |  2 )
    -0.1470905524127554,0.1149203802569079,0.0000000000000000,},
    {//(  4 |  3 )
    -0.0039797367270373,},
    {//(  4 |  4 )
    0.7600318835666090,},

};

TEST_CASE("Testing the Kinetic class"){

    auto molecule=make_molecule();
    auto bs=molecule.get_basis("sto-3g_cart");
    Integrals::LibIntKinetic<> mod;

    auto T = mod.run(molecule, {bs, bs});
      //    nwx_libint::Kinetic libints(0,molecule,bs,bs);
//    Integrals::TwoCenterIntegral *Ints = &libints;
//
//    size_t counter=0;
//    for(size_t si=0;si<5;++si)
//    {
//        const size_t ni = bs.shellsize(si);
//        for(size_t sj=0;sj<=si;++sj)
//        {
//            const size_t nj = bs.shellsize(sj);
//            std::vector<const double*> buf_vec=Ints->calculate(si,sj);
//            auto buffer = buf_vec[0];
//            if(buffer==nullptr)continue;
//            ptr_wrapper wrapped_buffer={buffer,ni*nj};
//            compare_integrals(wrapped_buffer, corr[counter]);
//            ++counter;
//        }
//    }
}

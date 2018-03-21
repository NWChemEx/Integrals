#include <Integrals/nwx_libint/LibInt2C.hpp>

#include "TestCommon.hpp"

//Computes the nuclear-electron energy integrals for water in STO-3G
std::vector<std::vector<double>> corr={
    {//(  0 |  0 )
    -61.5805952694322016,},
    {//(  1 |  0 )
    -7.4108218563311627,},
    {//(  1 |  1 )
    -10.0090711420700238,},
    {//(  2 |  0 )
    0.0000000000000000,-0.0144738837457361,0.0000000000000000,},
    {//(  2 |  1 )
    0.0000000000000000,-0.1768908347336429,0.0000000000000000,},
    {//(  2 |  2 )
    -9.9875499350885519,0.0000000000000000,0.0000000000000000,
     0.0000000000000000,-9.9440433416987553,0.0000000000000000,
     0.0000000000000000,0.0000000000000000,-9.8758759950909436,},
    {//(  3 |  0 )
    -1.2316855721424860,},
    {//(  3 |  1 )
    -2.9772268535781317,},
    {//(  3 |  2 )
    -1.8222369134761314,-1.4717933387129600,0.0000000000000000,},
    {//(  3 |  3 )
    -5.3002032522950184,},
    {//(  4 |  0 )
    -1.2316855721424860,},
    {//(  4 |  1 )
    -2.9772268535781312,},
    {//(  4 |  2 )
    1.8222369134761311,-1.4717933387129600,0.0000000000000000,},
    {//(  4 |  3 )
    -1.0671710804724346,},
    {//(  4 |  4 )
    -5.3002032522950167,},
};

TEST_CASE("Testing NuclearElectron class"){

    auto molecule=make_molecule();
    auto bs=molecule.get_basis("sto-3g_cart");
    nwx_libint::NuclearElectron libints(0,molecule,bs,bs);
    Integrals::TwoCenterIntegral *Ints = &libints;

    size_t counter=0;
    for(size_t si=0;si<5;++si)
    {
        const size_t ni = bs.shellsize(si);
        for(size_t sj=0;sj<=si;++sj)
        {
            const size_t nj = bs.shellsize(sj);
            std::vector<const double*> buf_vec=Ints->calculate(si,sj);
            auto buffer = buf_vec[0];
            if(buffer==nullptr)continue;
            ptr_wrapper wrapped_buffer={buffer,ni*nj};
            compare_integrals(wrapped_buffer, corr[counter]);
            ++counter;
        }
    }
}

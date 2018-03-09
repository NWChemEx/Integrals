#include <Integrals/nwx_libint/LibInt2C.hpp>

#include "TestCommon.hpp"

//Computes the density fitting metric integrals for water in STO-3G
//Note: the integrals actually use STO-3G and not a fitting basis, but I fail
//to see how that really matters for a unit test...
std::vector<std::vector<double>> corr={
    { //( 0 | 0 )
    1.0464370899978459, },
    { //( 1 | 0 )
    3.4291996305312606, },
    { //( 1 | 1 )
    26.4352252164276713, },
    { //( 2 | 0 )
    0.0000000000000000, 0.0000000000000000, 0.0000000000000000, },
    { //( 2 | 1 )
    0.0000000000000000, -0.0000000000000002, 0.0000000000000000, },
    { //( 2 | 2 )
      5.7847978365504300, 0.0000000000000000, 0.0000000000000000,
      0.0000000000000000, 5.7847978365504300, 0.0000000000000000,
      0.0000000000000000, 0.0000000000000000, 5.7847978365504300, },
    { //( 3 | 0 )
    2.6052624057150808, },
    { //( 3 | 1 )
    25.3420821293274017, },
    { //( 3 | 2 )
    4.2141100538676941, 3.2924421173969129, 0.0000000000000000, },
    { //( 3 | 3 )
    39.9325707858561643, },
    { //( 4 | 0 )
    2.6052624057150808, },
    { //( 4 | 1 )
    25.3420821293274017, },
    { //( 4 | 2 )
    -4.2141100538676941, 3.2924421173969129, 0.0000000000000000, },
    { //( 4 | 3 )
    26.6712894368539963, },
    { //( 4 | 4 )
    39.9325707858561643, },

};

TEST_CASE("Testing the Metric class"){

    auto atoms=make_atoms();
    auto bs=get_basis("PRIMARY",atoms);
    nwx_libint::Metric libints(0,atoms,bs,bs);
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

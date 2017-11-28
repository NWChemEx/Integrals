#include <IntegralsEx/nwx_libint/LibInt2C.hpp>
#define CATCH_CONFIG_MAIN
#include "TestCommon.hpp"

//Computes the density fitting metric integrals for water in STO-3G
//Note: the integrals actually use STO-3G and not a fitting basis, but I fail
//to see how that really matters for a unit test...
std::vector<std::vector<double>> corr={
    { //( 0 | 0 )
    1.0464371073132783, },
    { //( 1 | 0 )
    3.4291997049832545, },
    { //( 1 | 1 )
    26.4352259268832768, },
    { //( 2 | 0 )
    0.0000000000000000, 0.0000000000000000, 0.0000000000000000, },
    { //( 2 | 1 )
    0.0000000000000000, -0.0000000000000002, 0.0000000000000000, },
    { //( 2 | 2 )
      5.7847979494771771, 0.0000000000000000, 0.0000000000000000,
      0.0000000000000000, 5.7847979494771771, 0.0000000000000000,
      0.0000000000000000, 0.0000000000000000, 5.7847979494771771, },
    { //( 3 | 0 )
    2.6052624154025357, },
    { //( 3 | 1 )
    25.3420823544301896, },
    { //( 3 | 2 )
    4.2141100758044443, 3.2924421345358774, 0.0000000000000000, },
    { //( 3 | 3 )
    39.9325704220624615, },
    { //( 4 | 0 )
    2.6052624154025357, },
    { //( 4 | 1 )
    25.3420823544301896, },
    { //( 4 | 2 )
    -4.2141100758044443, 3.2924421345358774, 0.0000000000000000, },
    { //( 4 | 3 )
    26.6712891938732142, },
    { //( 4 | 4 )
    39.9325704220624615, },

};

TEST_CASE("Testing the Metric class"){

    auto atoms=make_atoms();
    auto bs=get_basis("PRIMARY",atoms);
    nwx_libint::Metric libints(0,atoms,bs,bs);
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

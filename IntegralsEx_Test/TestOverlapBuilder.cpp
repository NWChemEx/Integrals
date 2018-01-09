#include <IntegralsEx/TwoCTensorBuilder.hpp>
#define CATCH_CONFIG_MAIN
#include "TestCommon.hpp"

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

TEST_CASE("Testing OverlapTensorBuilder"){
    
    auto atoms=make_atoms();
    auto bs=get_basis("PRIMARY",atoms);

    std::vector<LibChemist::BasisSet> basissets({bs,bs});    

    IntegralsEx::TwoCTensorBuilder<nwx_libint::Overlap> overbuild;
    auto overtensor = overbuild.compute(atoms,basissets);

    REQUIRE(overtensor[0](0,0) == Approx(corr[0][0]));
    REQUIRE(overtensor[0](1,0) == Approx(corr[1][0]));
    
    std::vector<LibChemist::BasisSet> badsets({bs,bs,bs});  
    REQUIRE_THROWS_AS(overbuild.compute(atoms,badsets), std::length_error);  
}

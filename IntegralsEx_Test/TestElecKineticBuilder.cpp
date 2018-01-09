#include <IntegralsEx/TwoCTensorBuilder.hpp>
#define CATCH_CONFIG_MAIN
#include "TestCommon.hpp"

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

TEST_CASE("Testing ElecKineticTensorBuilder"){
    
    auto atoms=make_atoms();
    auto bs=get_basis("PRIMARY",atoms);

    std::vector<LibChemist::BasisSet> basissets({bs,bs});    

    IntegralsEx::TwoCTensorBuilder<nwx_libint::Kinetic> kinetic_build;
    auto kinetic_tensor = kinetic_build.compute(atoms,basissets);

    REQUIRE(kinetic_tensor[0](0,0) == Approx(corr[0][0]));
    REQUIRE(kinetic_tensor[0](1,0) == Approx(corr[1][0]));
    
    std::vector<LibChemist::BasisSet> badsets({bs,bs,bs});  
    REQUIRE_THROWS_AS(kinetic_build.compute(atoms,badsets), std::length_error);  
}

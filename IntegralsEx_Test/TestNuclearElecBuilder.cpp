#include <IntegralsEx/TwoCTensorBuilder.hpp>
#define CATCH_CONFIG_MAIN
#include "TestCommon.hpp"

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

TEST_CASE("Testing NuclearElecTensorBuilder"){
    
    auto atoms=make_atoms();
    auto bs=get_basis("PRIMARY",atoms);

    std::vector<LibChemist::BasisSet> basissets({bs,bs});    

    IntegralsEx::TwoCTensorBuilder<nwx_libint::NuclearElectron> nuclelec_build;
    auto nuclelec_tensor = nuclelec_build.compute(atoms,basissets);

    REQUIRE(nuclelec_tensor[0](0,0) == Approx(corr[0][0]));
    REQUIRE(nuclelec_tensor[0](1,0) == Approx(corr[1][0]));
    
    std::vector<LibChemist::BasisSet> badsets({bs,bs,bs});  
    REQUIRE_THROWS_AS(nuclelec_build.compute(atoms,badsets), std::length_error);  
}

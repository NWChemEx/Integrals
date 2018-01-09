#include <IntegralsEx/ThreeCTensorBuilder.hpp>
#define CATCH_CONFIG_MAIN
#include "TestCommon.hpp"
#include "H2O_STO3G_DF.hpp"

TEST_CASE("Testing DFTensorBuilder"){
    
    auto atoms=make_atoms();
    auto bs=get_basis("PRIMARY",atoms);

    std::vector<LibChemist::BasisSet> basissets({bs,bs,bs});    

    IntegralsEx::ThreeCTensorBuilder<nwx_libint::DF3C2E> df_build;
    auto df_tensor = df_build.compute(atoms,basissets);

    REQUIRE(df_tensor[0](0,0,0) == Approx(corr[0][0]));
    REQUIRE(df_tensor[0](0,0,1) == Approx(corr[1][0]));
    
    std::vector<LibChemist::BasisSet> badsets({bs,bs,bs,bs});  
    REQUIRE_THROWS_AS(df_build.compute(atoms,badsets), std::length_error);  
}

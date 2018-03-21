#include <Integrals/ThreeCTensorBuilder.hpp>

#include "TestCommon.hpp"
#include "H2O_STO3G_DF.hpp"

TEST_CASE("Testing DFTensorBuilder"){
    
    auto molecule=make_molecule();
    auto bs=molecule.get_basis("sto-3g_cart");

    std::vector<LibChemist::BasisSet> basissets({bs,bs,bs});    

    Integrals::ThreeCTensorBuilder<nwx_libint::DF3C2E> df_build;
    auto df_tensor = df_build.compute(molecule,basissets);

    size_t counter = 0;
    for (size_t i = 0; i < df_tensor[0].dimension(0); i++)
        for (size_t j = 0; j < df_tensor[0].dimension(1); j++) 
            for (size_t k = 0; k < df_tensor[0].dimension(2); k++) {
                REQUIRE(df_tensor[0](i,j,k) == Approx(corr_big[counter]).epsilon(eps).margin(marg));
                counter++;
        }
    
    std::vector<LibChemist::BasisSet> badsets({bs,bs,bs,bs});  
    REQUIRE_THROWS_AS(df_build.compute(molecule,badsets), std::length_error);  
}

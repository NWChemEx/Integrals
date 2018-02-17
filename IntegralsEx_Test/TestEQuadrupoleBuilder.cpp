#include <IntegralsEx/TwoCTensorBuilder.hpp>
#include "TestCommon.hpp"
#include "H2O_STO3G_Multipole.hpp"

TEST_CASE("Testing EQuadrupoleTensorBuilder"){
    
    auto atoms=make_atoms();
    auto bs=get_basis("PRIMARY",atoms);

    std::vector<LibChemist::BasisSet> basissets({bs,bs});    

    IntegralsEx::TwoCTensorBuilder<nwx_libint::EQuadrupole> quadrupole_build;
    auto quadrupole_tensor = quadrupole_build.compute(atoms,basissets);

    size_t counter = 0;
    for (size_t i = 0; i < quadrupole_tensor[0].dimension(0); i++)
        for (size_t j = 0; j < quadrupole_tensor[0].dimension(1); j++) {
            for (size_t ncomp = 0; ncomp < 10; ncomp++)
                REQUIRE(quadrupole_tensor[ncomp](i,j) == Approx(corr[ncomp][counter]).epsilon(eps).margin(marg));
            counter++;
        }
    
    std::vector<LibChemist::BasisSet> badsets({bs,bs,bs});  
    REQUIRE_THROWS_AS(quadrupole_build.compute(atoms,badsets), std::length_error);  
}

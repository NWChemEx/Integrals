#include <Integrals/TwoCTensorBuilder.hpp>
#include "TestCommon.hpp"
#include "H2O_STO3G_Multipole.hpp"

TEST_CASE("Testing EDipoleTensorBuilder"){
    
    auto molecule=make_molecule();
    auto bs=molecule.get_basis("sto-3g");

    std::vector<LibChemist::BasisSet> basissets({bs,bs});    

    Integrals::TwoCTensorBuilder<nwx_libint::EDipole> dipole_build;
    auto dipole_tensor = dipole_build.compute(molecule,basissets);

    size_t counter = 0;
    for (size_t i = 0; i < dipole_tensor[0].dimension(0); i++)
        for (size_t j = 0; j < dipole_tensor[0].dimension(1); j++) {
            for (size_t ncomp = 0; ncomp < 4; ncomp++)
                REQUIRE(dipole_tensor[ncomp](i,j) == Approx(corr[ncomp][counter]).epsilon(eps).margin(marg));
            counter++;
        }
    
    std::vector<LibChemist::BasisSet> badsets({bs,bs,bs});  
    REQUIRE_THROWS_AS(dipole_build.compute(molecule,badsets), std::length_error);  
}

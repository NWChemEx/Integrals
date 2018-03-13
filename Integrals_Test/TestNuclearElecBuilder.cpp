#include <Integrals/TwoCTensorBuilder.hpp>

#include "TestCommon.hpp"

std::vector<double> corr={
-61.5805952694322016,
-7.4108218563311627,
0.0000000000000000,
-0.0144738837457361,
0.0000000000000000,
-1.2316855721424858,
-1.2316855721424858,
-7.4108218563311627,
-10.0090711420700238,
0.0000000000000000,
-0.1768908347336429,
0.0000000000000000,
-2.9772268535781317,
-2.9772268535781312,
0.0000000000000000,
0.0000000000000000,
-9.9875499350885519,
0.0000000000000000,
0.0000000000000000,
-1.8222369134761285,
1.8222369134761285,
-0.0144738837457361,
-0.1768908347336429,
0.0000000000000000,
-9.9440433416987553,
0.0000000000000000,
-1.4717933387129591,
-1.4717933387129591,
0.0000000000000000,
0.0000000000000000,
0.0000000000000000,
0.0000000000000000,
-9.8758759950909436,
0.0000000000000000,
0.0000000000000000,
-1.2316855721424860,
-2.9772268535781317,
-1.8222369134761287,
-1.4717933387129594,
0.0000000000000000,
-5.3002032522950184,
-1.0671710804724346,
-1.2316855721424860,
-2.9772268535781312,
1.8222369134761287,
-1.4717933387129594,
0.0000000000000000,
-1.0671710804724346,
-5.3002032522950167
};

TEST_CASE("Testing NuclearElecTensorBuilder"){
    
    auto molecule=make_molecule();
    auto bs=molecule.get_basis("sto-3g");

    std::vector<LibChemist::BasisSet> basissets({bs,bs});    

    Integrals::TwoCTensorBuilder<nwx_libint::NuclearElectron> nuclelec_build;
    auto nuclelec_tensor = nuclelec_build.compute(molecule,basissets);

    size_t counter = 0;
    for (size_t i = 0; i < nuclelec_tensor[0].dimension(0); i++)
        for (size_t j = 0; j < nuclelec_tensor[0].dimension(1); j++) {
            REQUIRE(nuclelec_tensor[0](i,j) == Approx(corr[counter]).epsilon(eps).margin(marg));
            counter++;
        }

    std::vector<LibChemist::BasisSet> badsets({bs,bs,bs});  
    REQUIRE_THROWS_AS(nuclelec_build.compute(molecule,badsets), std::length_error);  
}

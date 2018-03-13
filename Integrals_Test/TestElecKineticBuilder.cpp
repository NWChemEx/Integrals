#include <Integrals/TwoCTensorBuilder.hpp>

#include "TestCommon.hpp"

std::vector<double> corr={
29.0031999455395848,
-0.1680109393164922,
0.0000000000000000,
0.0000000000000001,
0.0000000000000000,
-0.0084163851854474,
-0.0084163851854474,
-0.1680109393164923,
0.8081279549303477,
0.0000000000000000,
-0.0000000000000000,
0.0000000000000000,
0.0705173385189986,
0.0705173385189986,
0.0000000000000000,
0.0000000000000000,
2.5287311981947642,
0.0000000000000000,
0.0000000000000000,
0.1470905524127554,
-0.1470905524127554,
0.0000000000000001,
-0.0000000000000000,
0.0000000000000000,
2.5287311981947642,
0.0000000000000000,
0.1149203802569079,
0.1149203802569079,
0.0000000000000000,
0.0000000000000000,
0.0000000000000000,
0.0000000000000000,
2.5287311981947642,
0.0000000000000000,
0.0000000000000000,
-0.0084163851854474,
0.0705173385189986,
0.1470905524127554,
0.1149203802569079,
0.0000000000000000,
0.7600318835666090,
-0.0039797367270373,
-0.0084163851854474,
0.0705173385189986,
-0.1470905524127554,
0.1149203802569079,
0.0000000000000000,
-0.0039797367270373,
0.7600318835666090
};

TEST_CASE("Testing ElecKineticTensorBuilder"){
    
    auto molecule=make_molecule();
    auto bs=molecule.get_basis("sto-3g");

    std::vector<LibChemist::BasisSet> basissets({bs,bs});    

    Integrals::TwoCTensorBuilder<nwx_libint::Kinetic> kinetic_build;
    auto kinetic_tensor = kinetic_build.compute(molecule,basissets);

    size_t counter = 0;
    for (size_t i = 0; i < kinetic_tensor[0].dimension(0); i++)
        for (size_t j = 0; j < kinetic_tensor[0].dimension(1); j++) {
            REQUIRE(kinetic_tensor[0](i,j) == Approx(corr[counter]).epsilon(eps).margin(marg));
            counter++;
        }
    
    std::vector<LibChemist::BasisSet> badsets({bs,bs,bs});  
    REQUIRE_THROWS_AS(kinetic_build.compute(molecule,badsets), std::length_error);  
}

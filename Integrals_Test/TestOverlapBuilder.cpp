#include <Integrals/TwoCTensorBuilder.hpp>

#include "TestCommon.hpp"

std::vector<double> corr={
1.0000000000000004,
0.2367039365108476,
0.0000000000000000,
0.0000000000000000,
0.0000000000000000,
0.0384055905135490,
0.0384055905135490,
0.2367039365108476,
1.0000000000000002,
0.0000000000000000,
-0.0000000000000000,
0.0000000000000000,
0.3861387813310925,
0.3861387813310925,
0.0000000000000000,
0.0000000000000000,
1.0000000000000002,
0.0000000000000000,
0.0000000000000000,
0.2684376412681760,
-0.2684376412681760,
0.0000000000000000,
-0.0000000000000000,
0.0000000000000000,
1.0000000000000002,
0.0000000000000000,
0.2097276494226496,
0.2097276494226496,
0.0000000000000000,
0.0000000000000000,
0.0000000000000000,
0.0000000000000000,
1.0000000000000002,
0.0000000000000000,
0.0000000000000000,
0.0384055905135490,
0.3861387813310925,
0.2684376412681760,
0.2097276494226497,
0.0000000000000000,
1.0000000000000002,
0.1817608668218927,
0.0384055905135490,
0.3861387813310925,
-0.2684376412681760,
0.2097276494226497,
0.0000000000000000,
0.1817608668218927,
1.0000000000000002
};

TEST_CASE("Testing OverlapTensorBuilder"){
    
    auto molecule=make_molecule();
    auto bs=molecule.get_basis("sto-3g");

    std::vector<LibChemist::BasisSet> basissets({bs,bs});    

    Integrals::TwoCTensorBuilder<nwx_libint::Overlap> over_build;
    auto over_tensor = over_build.compute(molecule,basissets);

    size_t counter = 0;
    for (size_t i = 0; i < over_tensor[0].dimension(0); i++)
        for (size_t j = 0; j < over_tensor[0].dimension(1); j++) {
            REQUIRE(over_tensor[0](i,j) == Approx(corr[counter]).epsilon(eps).margin(marg));
            counter++;
        }
    
    std::vector<LibChemist::BasisSet> badsets({bs,bs,bs});  
    REQUIRE_THROWS_AS(over_build.compute(molecule,badsets), std::length_error);  
}
#include <Integrals/LibintIntegral.hpp>
#include "TestCommon.hpp"
#include "../Integrals/LibintIntegral.hpp"
#include "H2O_STO3G_Multipole.hpp" // holds the correct values

using namespace Integrals::Libint;

//Computes the quadrupole integrals for water in STO-3G

TEST_CASE("Testing Libint's Electric Quadrupole Integrals class"){
    auto [molecule, bs] = make_molecule();
    EQuadrupole EQuadBuilder;
    auto T = EQuadBuilder.run(molecule, {bs, bs});
    compare_integrals(T, corr);
}

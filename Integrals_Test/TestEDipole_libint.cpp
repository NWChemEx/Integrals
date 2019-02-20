#include <Integrals/LibintIntegral.hpp>
#include "TestCommon.hpp"
#include "H2O_STO3G_Multipole.hpp" // holds the correct values

using namespace Integrals::Libint;

//Computes the dipole integrals for water in STO-3G
TEST_CASE("Testing Libint's Electric Dipole Integrals class"){
    auto [molecule, bs] = make_molecule();
    EDipole EDipBuilder;
    auto T = EDipBuilder.run(molecule, {bs, bs});
    compare_integrals(T, corr);
}

#include <Integrals/LibintIntegral.hpp>
#include "TestCommon.hpp"
#include "../Integrals/LibintIntegral.hpp"
#include "H2O_STO3G_Multipole.hpp" // holds the correct values

using namespace Integrals::Libint;

//Computes the octopole integrals for water in STO-3G
TEST_CASE("Testing Libint's Electric Octopole Integrals class"){
    auto [molecule, bs] = make_molecule();
    EOctopole EOctBuilder;
    auto T = EOctBuilder.run(molecule, {bs, bs});
    compare_integrals(T, corr);
    //print_integrals(T);
}
#include <Integrals/LibintIntegral.hpp>
#include "TestCommon.hpp"
#include "H2O_STO3G_ERI.hpp"

using namespace Integrals::Libint;

//Computes the ERI integrals for water in STO-3G
TEST_CASE("Testing Libint's ERI"){
    auto [molecule, bs] = make_molecule();
    ERI IBuilder;
    auto Ints = IBuilder.run(molecule, {bs, bs, bs, bs});
    compare_integrals(Ints, corr);
}

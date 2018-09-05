#include <Integrals/LibintIntegral.hpp>
#include "TestCommon.hpp"
#include "H2O_STO3G_DF.hpp"

using namespace Integrals::Libint;

//Computes the three-center, two-electron integrals for water in sto-3g
//Note: as with the metric tensor using sto-3g instead of a proper fitting basis
//should be fine for testing purposes.

TEST_CASE("Testing DF3C2E"){
    auto [molecule, bs] = make_molecule();
    DFERI IBuilder;
    auto Ints = IBuilder.run(molecule, {bs, bs, bs});
    compare_integrals(Ints, corr);
}

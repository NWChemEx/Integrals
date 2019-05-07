#include <integrals/integralsmm.hpp>
#include <property_types/aointegral.hpp>
#include "TestCommon.hpp"
#include "H2O_STO3G_DF.hpp"

using namespace Integrals::Libint;

//Computes the three-center, two-electron integrals for water in sto-3g
//Note: as with the metric tensor using sto-3g instead of a proper fitting basis
//should be fine for testing purposes.

TEST_CASE("Testing DF3C2E"){
    using integral_type = property_types::AOIntegral<3, double>;
    SDE::ModuleManager mm;
    load_modules(mm);
    auto [molecule, bs] = make_molecule();
    std::array<LibChemist::AOBasisSet, 3> bases = {bs, bs, bs};
    auto [Ints] = mm.at("ERI3").run_as<integral_type>(molecule, bases, std::size_t{0});
    compare_integrals(Ints,corr);
}

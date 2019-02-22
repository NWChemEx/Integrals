#include <Integrals/LibintIntegral.hpp>
#include "LibChemist/AOBasisSet.hpp"
#include <LibChemist/Defaults/PropertyTypes.hpp>
#include <SDE/ModuleManager.hpp>
#include "TestCommon.hpp"
#include "H2O_STO3G_DF.hpp"

using namespace Integrals::Libint;

//Computes the three-center, two-electron integrals for water in sto-3g
//Note: as with the metric tensor using sto-3g instead of a proper fitting basis
//should be fine for testing purposes.

TEST_CASE("Testing DF3C2E"){
    SDE::ModuleManager mm;
    mm.add_module("Integral", std::make_shared<DFERI>());
    mm.set_default<LibChemist::AOIntegral>("Integral");
    auto [molecule, bs] = make_molecule();
    std::array<LibChemist::AOBasisSet, 3> basis_array = {bs, bs, bs};
    auto Ints = mm.at("Integral").run_as<LibChemist::AOIntegral>(molecule, basis_array);
/*  DFERI IBuilder;
    auto Ints = IBuilder.run(molecule, {bs, bs, bs});*/
    compare_integrals(Ints, corr);
}

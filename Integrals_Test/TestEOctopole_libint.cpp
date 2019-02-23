#include <Integrals/LibintIntegral.hpp>
#include <SDE/ModuleManager.hpp>
#include "TestCommon.hpp"
#include "H2O_STO3G_Multipole.hpp" // holds the correct values

using namespace Integrals::Libint;

//Computes the octopole integrals for water in STO-3G
TEST_CASE("Testing Libint's Electric Octopole Integrals class"){
    using integral_type = LibChemist::AOIntegral<2, double>;
    SDE::ModuleManager mm;
    mm.add_module("Integral", std::make_shared<EOctopole>());
    mm.set_default<integral_type>("Integral");
    auto [molecule, bs] = make_molecule();
    std::array<LibChemist::AOBasisSet, 2> bases = {bs, bs};
    auto [Ints] = mm.at("Integral").run_as<integral_type>(molecule, bases, std::size_t{0});
    compare_integrals(Ints, corr);
}
#include <Integrals/IntegralsMM.hpp>
#include <LibChemist/Defaults/PropertyTypes.hpp>
#include "TestCommon.hpp"
#include "H2O_STO3G_ERI.hpp"

using namespace Integrals::Libint;

//Computes the ERI integrals for water in STO-3G
TEST_CASE("Testing Libint's ERI"){
    using integral_type = LibChemist::AOIntegral<4, double>;
    SDE::ModuleManager mm;
    load_modules(mm);
    auto [molecule, bs] = make_molecule();
    std::array<LibChemist::AOBasisSet, 4> bases = {bs, bs, bs, bs};
    auto [Ints] = mm.at("ERI4").run_as<integral_type>(molecule, bases, std::size_t{0});
    compare_integrals(Ints, corr);
}

#include "H2O_STO3G_ERI.hpp"
#include "test_common.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/aointegral.hpp>

using namespace integrals::libint;

// Computes the ERI integrals for water in STO-3G
TEST_CASE("Testing Libint's ERI") {
    using integral_type = property_types::AOIntegral<4, double>;
    sde::ModuleManager mm;
    load_modules(mm);
    auto[molecule, bs]                          = make_molecule();
    std::array<libchemist::AOBasisSet, 4> bases = {bs, bs, bs, bs};
    auto[Ints] =
      mm.at("ERI4").run_as<integral_type>(molecule, bases, std::size_t{0});
    compare_integrals(Ints, corr);
}

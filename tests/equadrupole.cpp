#include "H2O_STO3G_Multipole.hpp" // holds the correct values
#include "test_common.hpp"
#include <integrals/integrals_mm.hpp>
#include <property_types/aointegral.hpp>

using namespace integrals::libint;

// Computes the quadrupole integrals for water in STO-3G

TEST_CASE("Testing Libint's Electric Quadrupole Integrals class") {
    using integral_type = property_types::AOIntegral<2, double>;
    sde::ModuleManager mm;
    load_modules(mm);
    auto[molecule, bs]                          = make_molecule();
    std::array<libchemist::AOBasisSet, 2> bases = {bs, bs};
    auto[Ints]                                  = mm.at("EQuadrupole")
                   .run_as<integral_type>(molecule, bases, std::size_t{0});
    compare_integrals(Ints, corr);
}

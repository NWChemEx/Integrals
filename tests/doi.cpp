#include "TestCommon.hpp"
#include "h2o_sto3g_doi.hpp"
#include <integrals/integrals_mm.hpp>
#include <property_types/aointegral.hpp>

TEST_CASE("DifferentialOverlapIntegrals") {
    using integral_type = property_types::AOIntegral<4, double>;
    sde::ModuleManager mm;
    integrals::libint::load_modules(mm);
    auto[molecule, bs]                          = make_molecule();
    std::array<libchemist::AOBasisSet, 4> bases = {bs, bs, bs, bs};
    auto[Ints] =
      mm.at("DOI").run_as<integral_type>(molecule, bases, std::size_t{0});
    compare_integrals(Ints, corr);
}

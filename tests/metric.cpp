#include "test_common.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/aointegral.hpp>

using namespace integrals::libint;

// Computes the density fitting metric integrals for water in STO-3G
// Note: the integrals actually use STO-3G and not a fitting basis, but I fail
// to see how that really matters for a unit test...
static BlockTensor corr{
  {{
     0,
     0,
   },
   {
     1.0464370899978459,  3.4291996305312606,  0.0000000000000000,
     0.0000000000000000,  0.0000000000000000,  2.6052624057150817,
     2.6052624057150817,  3.4291996305312606,  26.4352252164276713,
     -0.0000000000000002, 0.0000000000000000,  0.0000000000000000,
     25.3420821293274088, 25.3420821293274088, 0.0000000000000000,
     -0.0000000000000002, 5.7847978365504318,  0.0000000000000000,
     0.0000000000000000,  3.2924421173969143,  3.2924421173969143,
     0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
     5.7847978365504300,  0.0000000000000000,  0.0000000000000000,
     0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
     0.0000000000000000,  0.0000000000000000,  5.7847978365504318,
     4.2141100538676941,  -4.2141100538676941, 2.6052624057150817,
     25.3420821293274088, 3.2924421173969143,  0.0000000000000000,
     4.2141100538676941,  39.9325707858561643, 26.6712894368540034,
     2.6052624057150817,  25.3420821293274088, 3.2924421173969143,
     0.0000000000000000,  -4.2141100538676941, 26.6712894368540034,
     39.9325707858561643,
   }}};

TEST_CASE("Testing the Metric class") {
    using integral_type = property_types::AOIntegral<2, double>;
    sde::ModuleManager mm;
    load_modules(mm);
    auto[molecule, bs]                          = make_molecule();
    std::array<libchemist::AOBasisSet, 2> bases = {bs, bs};
    auto[Ints] =
      mm.at("ERI2").run_as<integral_type>(molecule, bases, std::size_t{0});
    compare_integrals(Ints, corr);
}

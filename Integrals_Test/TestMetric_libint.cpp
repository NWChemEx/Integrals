#include <Integrals/LibintIntegral.hpp>
#include <SDE/ModuleManager.hpp>
#include "TestCommon.hpp"

using namespace Integrals::Libint;

//Computes the density fitting metric integrals for water in STO-3G
//Note: the integrals actually use STO-3G and not a fitting basis, but I fail
//to see how that really matters for a unit test...
static BlockTensor corr{
{{0 , 0}, {1.0464370899978459,}},
{{0 , 1}, {3.4291996305312606,}},
{{0 , 2}, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
{{0 , 3}, {2.6052624057150817,}},
{{0 , 4}, {2.6052624057150817,}},
{{1 , 0}, {3.4291996305312606,}},
{{1 , 1}, {26.4352252164276713,}},
{{1 , 2}, {-0.0000000000000002,0.0000000000000000,0.0000000000000000,}},
{{1 , 3}, {25.3420821293274088,}},
{{1 , 4}, {25.3420821293274088,}},
{{2 , 0}, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
{{2 , 1}, {-0.0000000000000002,0.0000000000000000,0.0000000000000000,}},
{{2 , 2}, {5.7847978365504318,0.0000000000000000,0.0000000000000000,
           0.0000000000000000,5.7847978365504300,0.0000000000000000,
           0.0000000000000000,0.0000000000000000,5.7847978365504318,}},
{{2 , 3}, {3.2924421173969143,0.0000000000000000,4.2141100538676941,}},
{{2 , 4}, {3.2924421173969143,0.0000000000000000,-4.2141100538676941,}},
{{3 , 0}, {2.6052624057150817,}},
{{3 , 1}, {25.3420821293274088,}},
{{3 , 2}, {3.2924421173969143,0.0000000000000000,4.2141100538676941,}},
{{3 , 3}, {39.9325707858561643,}},
{{3 , 4}, {26.6712894368540034,}},
{{4 , 0}, {2.6052624057150817,}},
{{4 , 1}, {25.3420821293274088,}},
{{4 , 2}, {3.2924421173969143,0.0000000000000000,-4.2141100538676941,}},
{{4 , 3}, {26.6712894368540034,}},
{{4 , 4}, {39.9325707858561643,}}};

TEST_CASE("Testing the Metric class"){
    using integral_type = LibChemist::AOIntegral<2, double>;
    SDE::ModuleManager mm;
    mm.add_module("Integral", std::make_shared<Metric>());
    mm.set_default<integral_type>("Integral");
    auto [molecule, bs] = make_molecule();
    std::array<LibChemist::AOBasisSet, 2> bases = {bs, bs};
    auto [Ints] = mm.at("Integral").run_as<integral_type>(molecule, bases, std::size_t{0});
    compare_integrals(Ints, corr);
}

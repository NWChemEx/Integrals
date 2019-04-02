#include <Integrals/IntegralsMM.hpp>
#include <LibChemist/Defaults/PropertyTypes.hpp>
#include "TestCommon.hpp"

using namespace Integrals::Libint;

//Computes the overlap integrals for water in STO-3G
static BlockTensor corr{
{{0, 0, }, {
   1.0000000000000004,0.2367039365108476,0.0000000000000000,0.0000000000000000,0.0000000000000000,
   0.2367039365108476,1.0000000000000002,-0.0000000000000000,0.0000000000000000,0.0000000000000000,
   0.0000000000000000,-0.0000000000000000,1.0000000000000007,0.0000000000000000,0.0000000000000000,
   0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000002,0.0000000000000000,
   0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000007,}},
{{0, 1, }, {
   0.0384055905135491,0.3861387813310929,0.2097276494226498,0.0000000000000000,0.2684376412681763,}},
{{0, 2, }, {
   0.0384055905135491,0.3861387813310929,0.2097276494226498,0.0000000000000000,-0.2684376412681763,}},
{{1, 0, }, {
   0.0384055905135491,0.3861387813310928,0.2097276494226498,0.0000000000000000,0.2684376412681763,}},
{{1, 1, }, {
   1.0000000000000002,}},
{{1, 2, }, {
   0.1817608668218930,}},
{{2, 0, }, {
   0.0384055905135491,0.3861387813310928,0.2097276494226498,0.0000000000000000,-0.2684376412681763,}},
{{2, 1, }, {
   0.1817608668218930,}},
{{2, 2, }, {
   1.0000000000000002,}}
};

TEST_CASE("Testing LibIntOverlap class"){
    using integral_type = LibChemist::AOIntegral<2, double>;
    SDE::ModuleManager mm;
    load_modules(mm);
    mm.at("Overlap").change_input("Tile Size",size_t{1});
    auto [molecule, bs] = make_molecule();
    std::array<LibChemist::AOBasisSet, 2> bases = {bs, bs};
    auto [Ints] = mm.at("Overlap").run_as<integral_type>(molecule, bases, std::size_t{0});
    compare_integrals(Ints,corr);
}

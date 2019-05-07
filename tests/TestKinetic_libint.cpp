#include <integrals/integralsmm.hpp>
#include <integrals/libint_integral.hpp>
#include <property_types/aointegral.hpp>
#include "TestCommon.hpp"

using namespace Integrals::Libint;

//Computes the kinetic energy integrals for water in STO-3G
static BlockTensor corr{
{{0, 0, }, {
   29.0031999455395848,-0.1680109393164922,0.0000000000000000,0.0000000000000000,0.0000000000000000,
   -0.0084163851854474,-0.0084163851854474,-0.1680109393164923,0.8081279549303477,-0.0000000000000000,
   0.0000000000000000,0.0000000000000000,0.0705173385189988,0.0705173385189988,0.0000000000000000,
   -0.0000000000000000,2.5287311981947651,0.0000000000000000,0.0000000000000000,0.1149203802569082,
   0.1149203802569082,0.0000000000000000,0.0000000000000000,0.0000000000000000,2.5287311981947642,
   0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
   0.0000000000000000,0.0000000000000000,2.5287311981947651,0.1470905524127557,-0.1470905524127557,
   -0.0084163851854474,0.0705173385189988,0.1149203802569082,0.0000000000000000,0.1470905524127557,
   0.7600318835666090,-0.0039797367270372,-0.0084163851854474,0.0705173385189988,0.1149203802569082,
   0.0000000000000000,-0.1470905524127557,-0.0039797367270372,0.7600318835666090,}}};

TEST_CASE("Testing Libint's Kinetic Energy Integrals class"){
    auto[molecule, bs] = make_molecule();
    std::array<LibChemist::AOBasisSet, 2> bases = {bs, bs};
    using integral_type = property_types::AOIntegral<2, double>;
    SDE::ModuleManager mm;
    SECTION("Direct implementation") {
        load_modules(mm);
        auto[Ints] = mm.at("Kinetic").run_as<integral_type>(molecule, bases, std::size_t{0});
        compare_integrals(Ints,corr);
    }
    SECTION("Core implementation") {
        mm.add_module("KineticCore", std::make_shared<Kinetic>(detail_::implementation_type::core));
        auto[Ints] = mm.at("KineticCore").run_as<integral_type>(molecule, bases, std::size_t{0});
        compare_integrals(Ints,corr);
    }
}

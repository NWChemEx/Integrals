#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/yukawa.hpp>

TEST_CASE("Yukawa3C") {
    using integral_type = property_types::Yukawa3CIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    mm.at("Yukawa3").change_input("Tile size", std::size_t{1});
    auto [molecule, bs] = make_molecule();
    auto [I] = mm.at("Yukawa3").run_as<integral_type>(bs, bs, bs, std::size_t{0});

    std::cout << "Yukawa3" << std::endl;
    std::cout << I << std::endl;
    std::cout << std::endl;
}
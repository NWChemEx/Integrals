#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/stg.hpp>

TEST_CASE("STG2C") {
    using integral_type = property_types::STG2CIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    mm.at("STG2").change_input("Tile size", std::size_t{1});
    auto [molecule, bs] = make_molecule();
    auto [I] = mm.at("STG2").run_as<integral_type>(bs, bs, std::size_t{0});

    std::cout << "STG2" << std::endl;
    std::cout << I << std::endl;
    std::cout << std::endl;
}
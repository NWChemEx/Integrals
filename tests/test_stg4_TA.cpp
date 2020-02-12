#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/stg.hpp>

TEST_CASE("STG4C") {
    using integral_type = property_types::STG4CIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    mm.at("STG4").change_input("Tile size", std::vector<std::size_t>{1});
    auto [molecule, bs] = make_molecule();
    auto stg_exponent = 1.0;
    auto [I] = mm.at("STG4").run_as<integral_type>(bs, bs, bs, bs, std::size_t{0}, stg_exponent);

//    std::cout << "STG4" << std::endl;
//    std::cout << I << std::endl;
//    std::cout << std::endl;
}
#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/electron_repulsion.hpp>

TEST_CASE("ERI4C") {
    using integral_type = property_types::ERI4CIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    mm.at("ERI4").change_input("Tile size", std::size_t{1});
    auto [molecule, bs] = make_molecule();
    auto [I] = mm.at("ERI4").run_as<integral_type>(bs, bs, bs, bs, std::size_t{0});

    std::cout << "ERI4" << std::endl;
    std::cout << I << std::endl;
    std::cout << std::endl;
}

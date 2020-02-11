#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/kinetic.hpp>

TEST_CASE("Kinetic") {
    using integral_type = property_types::KineticIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    mm.at("Kinetic").change_input("Tile size", std::size_t{1});
    auto [molecule, bs] = make_molecule();
    auto [T] = mm.at("Kinetic").run_as<integral_type>(bs, bs, std::size_t{0});

    std::cout << "Kinetic" << std::endl;
    std::cout << T << std::endl;
    std::cout << std::endl;
}

#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/nuclear.hpp>

TEST_CASE("Nuclear") {
    using integral_type = property_types::NuclearIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    mm.at("Nuclear").change_input("Tile size", std::vector<std::size_t>{1});
    auto [molecule, bs] = make_molecule();
    auto [V] = mm.at("Nuclear").run_as<integral_type>(bs, bs, molecule, std::size_t{0});

//    std::cout << "Nuclear" << std::endl;
//    std::cout << V << std::endl;
//    std::cout << std::endl;
}

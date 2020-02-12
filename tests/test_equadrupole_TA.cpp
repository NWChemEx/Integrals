#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/emultipole.hpp>

TEST_CASE("Octopole") {
    using integral_type = property_types::EOctopoleIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    mm.at("EOctopole").change_input("Tile size", std::vector<std::size_t>{1});
    auto [molecule, bs] = make_molecule();
    auto origin = std::array<double, 3>{0,0,0};
    auto [X] = mm.at("EOctopole").run_as<integral_type>(bs, bs, std::size_t{0}, origin);

//    std::cout << "Octopole" << std::endl;
//    std::cout << X << std::endl;
//    std::cout << std::endl;
}
#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/emultipole.hpp>

TEST_CASE("Dipole") {
    using integral_type = property_types::EDipoleIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    mm.at("EDipole").change_input("Tile size", std::vector<std::size_t>{1});
    auto [molecule, bs] = make_molecule();
    auto origin = std::array<double, 3>{0,0,0};
    auto [X] = mm.at("EDipole").run_as<integral_type>(bs, bs, std::size_t{0}, origin);

//    std::cout << "Dipole" << std::endl;
//    std::cout << X << std::endl;
//    std::cout << std::endl;
}

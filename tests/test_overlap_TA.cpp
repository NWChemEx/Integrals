#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/overlap.hpp>

TEST_CASE("Overlap") {
    using integral_type = property_types::OverlapIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    mm.at("Overlap").change_input("Tile size", std::vector<std::size_t>{1});
    auto [molecule, bs] = make_molecule();
    auto [S] = mm.at("Overlap").run_as<integral_type>(bs, bs, std::size_t{0});

    std::cout << "Overlap" << std::endl;
    std::cout << S << std::endl;
    std::cout << std::endl;
}

#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/yukawa.hpp>
#include "H2O_STO3G_Yukawa[1]_3C.hpp"

TEST_CASE("Yukawa3C") {
    using integral_type = property_types::Yukawa3CIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto stg_exponent = 1.0;
    mm.at("Yukawa3").change_input("Tile size", std::vector<std::size_t>{1, 1, 4, 1});
    auto [I] = mm.at("Yukawa3").run_as<integral_type>(bs, bs, bs, std::size_t{0}, stg_exponent);

    compare_integrals(I, corr);
}
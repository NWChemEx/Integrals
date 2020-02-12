#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/yukawa.hpp>

TEST_CASE("Yukawa2C") {
    using integral_type = property_types::Yukawa2CIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto stg_exponent = 1.0;
    auto [I] = mm.at("Yukawa2").run_as<integral_type>(bs, bs, std::size_t{0}, stg_exponent);
}

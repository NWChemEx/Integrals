#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/yukawa.hpp>
#include "H2O_STO3G_Yukawa[1].hpp"

TEST_CASE("Yukawa4C") {
    using integral_type = property_types::Yukawa4CIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto stg_exponent = 1.0;
    mm.at("Yukawa4").change_input("Screening Threshold", 0.000001);
    auto [X] = mm.at("Yukawa4").run_as<integral_type>(bs, bs, bs, bs, std::size_t{0}, stg_exponent);

    compare_integrals(X, yukawa1ref, 0.0, 1e-12);
}

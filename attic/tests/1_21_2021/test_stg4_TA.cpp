#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/stg.hpp>
#include "H2O_STO3G_STG[1].hpp"

TEST_CASE("STG4C") {
    using integral_type = property_types::STG4CIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto stg_exponent = 1.0;
    mm.at("STG4").change_input("Screening Threshold", 0.000001);
    auto [X] = mm.at("STG4").run_as<integral_type>(bs, bs, bs, bs, std::size_t{0}, stg_exponent);

    REQUIRE(libchemist::ta_helpers::allclose(X, TensorType(X.world(), X.trange(), corr)));
}
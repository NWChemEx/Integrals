#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/electron_repulsion.hpp>
#include "H2O_STO3G_DF.hpp"

TEST_CASE("ERI3C") {
    using integral_type = property_types::ERI3CIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    mm.at("ERI3").change_input("Screening Threshold", 0.000001);
    auto [X] = mm.at("ERI3").run_as<integral_type>(bs, bs, bs, std::size_t{0});

    compare_integrals(X, corr);
}

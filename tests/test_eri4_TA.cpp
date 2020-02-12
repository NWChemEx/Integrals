#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/electron_repulsion.hpp>
#include "H2O_STO3G_ERI.hpp"

TEST_CASE("ERI4C") {
    using integral_type = property_types::ERI4CIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto [I] = mm.at("ERI4").run_as<integral_type>(bs, bs, bs, bs, std::size_t{0});

    compare_integrals(I, corr);
}

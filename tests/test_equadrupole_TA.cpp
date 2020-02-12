#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/emultipole.hpp>
#include "H2O_STO3G_Multipole.hpp"

TEST_CASE("Quadrupole") {
    using integral_type = property_types::EQuadrupoleIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto origin = std::array<double, 3>{0,0,0};
    auto [X] = mm.at("EQuadrupole").run_as<integral_type>(bs, bs, std::size_t{0}, origin);

    compare_integrals(X, corr);
}

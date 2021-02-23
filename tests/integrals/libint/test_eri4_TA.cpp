#include "../../test_common_TA.hpp"
#include "integrals/integrals.hpp"
#include "nwx_testing/H2O_STO3G_ERI.hpp"

TEST_CASE("ERI4C") {
    using integral_type = integrals::pt::eri4c<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto [X]            = mm.at("ERI4").run_as<integral_type>(bs, bs, bs, bs);

    REQUIRE(libchemist::ta_helpers::allclose(
      X, TensorType(X.world(), X.trange(), corr)));
}

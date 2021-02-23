#include "../../test_common_TA.hpp"
#include "integrals/integralsmm.hpp"
#include "integrals/property_types.hpp"
#include "nwx_testing/H2O_STO3G_DOI.hpp"
#include <libchemist/ta_helpers/ta_helpers.hpp>

TEST_CASE("DOI") {
    using integral_type = integrals::pt::doi<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto [X]            = mm.at("DOI").run_as<integral_type>(bs, bs);

    REQUIRE(libchemist::ta_helpers::allclose(
      X, TensorType(X.world(), X.trange(), corr)));
}
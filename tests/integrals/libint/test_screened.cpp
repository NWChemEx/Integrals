#include "../../tensors/H2O_STO3G_ERI_SCREENED.hpp"
#include "../../test_common_TA.hpp"
#include "integrals/integrals.hpp"

TEST_CASE("ERI4C CS") {
    using integral_type = integrals::pt::eri4c<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    mm.change_input("ERI4 CS", "Screening Threshold", 0.005);
    auto [molecule, bs] = make_molecule();
    auto [X]            = mm.at("ERI4 CS").run_as<integral_type>(bs, bs, bs, bs);

    REQUIRE(libchemist::ta_helpers::allclose(
      X, TensorType(X.world(), X.trange(), corr)));
}

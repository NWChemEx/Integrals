#include "../../tensors/H2O_STO3G_DF.hpp"
#include "../../test_common_TA.hpp"
#include "integrals/integrals.hpp"

TEST_CASE("ERI3C") {
    using integral_type = integrals::pt::eri3c<double>;
    const auto key      = "ERI3";

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    // mm.change_input(key, "Screening Threshold", 0.000001);
    auto [X] = mm.run_as<integral_type>(key, bs, bs, bs);
    TensorType X_corr(X.world(), X.trange(), corr);
    REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));
}

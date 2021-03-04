#include "../../tensors/H2O_STO3G_DF.hpp"
#include "../../test_common_TA.hpp"
#include "integrals/integrals.hpp"

TEST_CASE("ERI3C") {
    using integral_type = integrals::pt::eri3c<double>;
    const auto key1     = "ERI3";
    const auto key2     = "ERI3 CS";

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto [X] = mm.run_as<integral_type>(key1, bs, bs, bs);
    TensorType X_corr(X.world(), X.trange(), corr);
    REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));

    mm.change_input(key2, "Screening Threshold", 0.000001);
    auto [X2] = mm.run_as<integral_type>(key2, bs, bs, bs);
    REQUIRE(libchemist::ta_helpers::allclose(X2, X_corr));
}
